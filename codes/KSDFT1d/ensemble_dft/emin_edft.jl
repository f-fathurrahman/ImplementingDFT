function edft_calc_h_matrix!(Ham, psi, hmat)
    Npoints = size(psi, 1)
    Nstates = size(psi, 2)
    hx = Ham.grid.hx
    fill!(hmat, 0.0)
    # Apply T + V_ions to each orbital (column)
    H0 = Ham.Kmat + diagm( 0 => Ham.potentials.Ions)
    # XXX just use matmul?
    for jst in 1:Nstates
        Hpsi_j = H0 * psi[:,jst]
        for ist in 1:Nstates
            hmat[ist,jst] = dot(psi[:,ist], Hpsi_j)*hx
        end
    end
    return
end

function edft_calc_rhoe!(
    Ham::Hamiltonian1d,
    Focc_matrix::Matrix{Float64},
    psi::Matrix{Float64},
    rhoe::Matrix{Float64}
)
    ispin = 1
    Npoints = size(rhoe, 1)
    Nstates = size(psi, 2)
    Nelectrons = Ham.electrons.Nelectrons
    hx = Ham.grid.hx
    #
    contrib = zeros(Float64, Npoints)
    fill!(rhoe, 0.0)
    for ist in 1:Nstates
        for jst in ist:Nstates
            if ist == jst
                @. contrib[:] = Focc_matrix[ist,ist] * psi[:,ist] * psi[:,ist]
            else
                @. contrib[:] = 2.0 * Focc_matrix[ist,jst] * psi[:,ist] * psi[:,jst]
            end
            rhoe[:,ispin] .+= contrib[:]
        end
    end
    for i in 1:length(rhoe)
        if rhoe[i] < 0.0
            rhoe[i] = 0.0
        end
    end
    #ss = sum(rhoe)*hx
    ##println("ss = ", ss)
    #rhoe[:] *= Nelectrons/ss
    
    #ss = sum(rhoe)*hx
    #println("ss = ", ss)
    return
end

function edft_calc_Vhxc_mat!(Ham, psi, Vhxc, Vhxc_mat)
    Nstates = size(psi, 2)
    hx = Ham.grid.hx
    for ist in 1:Nstates, jst in ist:Nstates
        vij = dot( psi[:,ist] .* psi[:,jst], Vhxc) * hx
        Vhxc_mat[ist,jst] = vij
        if ist != jst
            Vhxc_mat[jst,ist] = vij
        end
    end
    return
end

function ortho_gram_schmidt!(psi, hx)
    Nstates = size(psi, 2)
    for ist in 1:Nstates
        nrm = sqrt(dot(psi[:,ist], psi[:,ist])*hx)
        if nrm > 1e-12
            psi[:,ist] /= nrm
        end
        for jst in (ist+1):Nstates
            prj = dot(psi[:,jst], psi[:,ist])*hx
            psi[:,jst] .-= prj * psi[:,ist]
        end
    end
    return
end


"""
Return S and s'(f) matrix in the current basis.
For degeneracy g, s(f) = -[ f ln(f/g) + (g-f) ln(1 - f/g) ].
"""
function calc_entropy(Focc_matrix)

    # HARDCODED
    ispin = 1
    degeneracy = 2

    f_diag, X = eigen(Hermitian(Focc_matrix))          # evecs are columns of unitary
    SMALL = 1e-12
    # Clamp eigenvalues to [0, g]
    f_diag = clamp.(f_diag, SMALL, degeneracy - SMALL)
    #
    S = -sum( @. f_diag * log(f_diag/degeneracy) + (degeneracy - f_diag) * log(1 - f_diag/degeneracy) )
    #
    sprim_diag = @. -log(f_diag / (degeneracy - f_diag))   # derivative d s / d f_i
    sprim = X * diagm(0 => sprim_diag) * X'
    return S, sprim
end


function calc_free_energy(Ham, Focc_matrix, hmat, Vhartree, rhoe)
    hx = Ham.grid.hx
    E_kext = tr(Focc_matrix * hmat)
    E_H = 0.5 * dot(rhoe, Vhartree) * hx
    
    epsxc = calc_epsxc_1d(Ham.xc_calc, rhoe[:,1])
    E_xc = dot(rhoe, epsxc)*hx

    kT = Ham.electrons.kT
    S, _ = calc_entropy(Focc_matrix)
    A = E_kext + E_H + E_xc - kT*S
    return A
end

function inner_loop!(
    Ham, Focc_matrix, psi, hmat;
    TOL = 1e-8,
    NmaxIter_inner = 2,
    verbose = true
)

    ebands = Ham.electrons.ebands
    Nspin = Ham.electrons.Nspin
    Npoints = Ham.grid.Npoints
    Vhartree = Ham.potentials.Hartree
    Vxc = Ham.potentials.XC
    Vtot = Ham.potentials.Total
    Nstates = Ham.electrons.Nstates
    Nelectrons = Ham.electrons.Nelectrons
    kT = Ham.electrons.kT

    if Nspin == 1
        degeneracy = 2
    else
        degeneracy = 1
    end
    ispin = 1  # HARDCODED

    Vhxc_mat = zeros(Float64, Nstates, Nstates)
    Hsub = zeros(Float64, Nstates, Nstates)
    rhoe_current = zeros(Float64, Npoints, Nspin)

    rhoe_1 = zeros(Float64, Npoints, Nspin)
    Vhartree1 = similar(Vhartree)
    Vxc1 = similar(Vxc)
    Vhxc_mat1 = zeros(Float64, Nstates, Nstates)
    Vtot1 = similar(Vtot)

    rhoe_β = zeros(Float64, Npoints, Nspin)
    Vhartree_β = similar(Vhartree)
    Vxc_β = similar(Vxc)

    rhoe_final = zeros(Float64, Npoints, Nspin)

    change = Inf # default?

    for iter_inner in 1:NmaxIter_inner

        edft_calc_rhoe!(Ham, Focc_matrix, psi, rhoe_current)

        ρ = reshape(rhoe_current, Npoints) # should be sum over Nspin
        Poisson_solve_sum!(Ham.grid, ρ, Vhartree)
        Vxc[:,ispin] = calc_Vxc_1d(Ham.xc_calc, rhoe_current)
        Vtot[:,ispin] = Vhartree[:] + Vxc[:,ispin]

        edft_calc_Vhxc_mat!(Ham, psi, Vtot, Vhxc_mat)
        Hsub[:,:] = hmat + Vhxc_mat

        ebands[:,ispin], X = eigen(Hermitian(Hsub))

        E_f = find_E_fermi( smear_fermi, ebands, Nelectrons, kT )
        # Stable Fermi‑Dirac using expit
        Focc_diag = @. degeneracy * expit(-(ebands[:,1] - E_f) / kT)
        Ham.electrons.Focc[:,ispin] = Focc_diag[:]

        Focc_matrix_tilde = X * diagm(0 => Focc_diag) * X'
        #Focc_matrix_tilde[:,:] = 0.5*(Focc_matrix_tilde + Focc_matrix_tilde')

        delta_f = Focc_matrix_tilde - Focc_matrix
        change = norm(delta_f)

        A0 = calc_free_energy(Ham, Focc_matrix, hmat, Vhartree, rhoe_current)
    
        edft_calc_rhoe!(Ham, Focc_matrix_tilde, psi, rhoe_1)
        ρ = reshape(rhoe_1, Npoints) # should be sum over Nspin
        Poisson_solve_sum!(Ham.grid, ρ, Vhartree1)
        Vxc1[:,ispin] = calc_Vxc_1d(Ham.xc_calc, rhoe_1)

        A1 = calc_free_energy(Ham, Focc_matrix_tilde, hmat, Vhartree1, rhoe_1)

        _, sprim0 = calc_entropy(Focc_matrix)
        dA0 = sum(delta_f .* (Hsub - kT*sprim0))

        _, sprim1 = calc_entropy(Focc_matrix_tilde)


        Vtot1[:,ispin] = Vhartree1[:] + Vxc1[:,ispin]

        edft_calc_Vhxc_mat!(Ham, psi, Vtot1, Vhxc_mat1)
        Hsub1 = hmat + Vhxc_mat1
        dA1 = sum(delta_f .* (Hsub1 - kT*sprim1))

        # Cubic interpolation for beta
        a = dA1 + dA0 - 2*(A1 - A0)
        b = 3*(A1 - A0) - 2*dA0 - dA1
        c = dA0
        d = A0

        β_opt = NaN
        A_opt = Inf

        if abs(a) > 1e-12
            disc = b*b - 3*a*c
            if disc >= 0
                sqrt_disc = sqrt(disc)
                β1 = (-b + sqrt_disc) / (3*a)
                β2 = (-b - sqrt_disc) / (3*a)
                for β in (β1, β2)
                    if 0 <= β <= 1
                        Focc_matrix_β = Focc_matrix + β * delta_f
                        edft_calc_rhoe!(Ham, Focc_matrix_β, psi, rhoe_β)
                        ρ = reshape(rhoe_β, Npoints) # should be sum over Nspin
                        Poisson_solve_sum!(Ham.grid, ρ, Vhartree_β)
                        Vxc_β[:,ispin] = calc_Vxc_1d(Ham.xc_calc, rhoe_β)
                        #
                        A_β = calc_free_energy(Ham, Focc_matrix_β, hmat, Vhartree_β, rhoe_β)
                        #
                        if A_β < A_opt
                            A_opt = A_β
                            β_opt = β
                        end
                    end
                end
            end
        end
        
        println("β_opt = $(β_opt) A_opt = $(A_opt)")

        if !isnan(β_opt)
            # Update f
            @views Focc_matrix[:,:] .+= β_opt * delta_f
            #
            verbose && println("      inner $(iter_inner): β = $(β_opt), ΔA = $(A_opt - A0)")
        else
            if A1 < A0
                @views Focc_matrix[:,:] = Focc_matrix_tilde[:,:]
                #
                verbose && println("      inner $(iter_inner): using β=1, ΔA = $(A1 - A0)")
            else
                verbose && println("      inner $(iter_inner): no improvement, β=0 kept")
            end
        end

        if change < TOL
            verbose && println("Inner loop converged")
            break
        end
    end

    edft_calc_rhoe!(Ham, Focc_matrix, psi, rhoe_final)

    return change, rhoe_final

end

function density_difference(rhoe_1, rhoe_2, hx)
    return sum(abs.(rhoe_1 - rhoe_2)) * hx
end


function emin_edft!(Ham)
    
    grid = Ham.grid
    Npoints = grid.Npoints
    Nstates = Ham.electrons.Nstates
    hx = grid.hx
    Nspin = Ham.electrons.Nspin

    Kmat = Ham.Kmat
    Vtot = Ham.potentials.Total
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree

    psi = zeros(Float64, Npoints, Nstates)
    initial_wf_mode = :gaussian
    if initial_wf_mode == :gaussian
        centers = range(grid.x[1], grid.x[end], length=Nstates)
        println("centers = ", collect(centers))
        width = (grid.x[end] - grid.x[1]) / Nstates
        println("width = ", width)
        for ist in 1:Nstates, ip in 1:Npoints
            psi[ip,ist] = exp(-((grid.x[ip] - centers[ist])/width)^2)
        end
    elseif mode == :random
        psi[:,:] = randn(Npoints, Nstates)
    else
        error("mode must be :gaussian or :random")
    end

    #ortho_sqrt!(psi)
    #psi[:,:] = psi[:,:]/sqrt(hx)
    
    ortho_gram_schmidt!(psi, hx)

    hmat = zeros(Float64, Nstates, Nstates)
    edft_calc_h_matrix!(Ham, psi, hmat)
    display(hmat); println()

    Nelectrons = Ham.electrons.Nelectrons
    kT = Ham.electrons.kT
    ebands = Ham.electrons.ebands
    ispin = 1
    
    ebands[:,1], X = eigen(Hermitian(hmat))
    E_f = find_E_fermi( smear_fermi, ebands, Nelectrons, kT )
    println("E_f = ", E_f)

    if Nspin == 1
        degeneracy = 2
    else
        degeneracy = 1
    end

    Focc_diag = @. degeneracy * expit(-(ebands[:,1] - E_f) / kT)
    println("Focc_diag = ", Focc_diag)

    Focc_matrix = X * diagm(0 => Focc_diag) * X'
    Focc_matrix[:,:] = 0.5*(Focc_matrix + Focc_matrix')

    # Begin outer loop

    rhoe = zeros(Float64, Npoints, Nspin)
    rhoe_old = zeros(Float64, Npoints, Nspin)

    edft_calc_rhoe!(Ham, Focc_matrix, psi, rhoe_old)
    #
    ρ = reshape(rhoe_old, Npoints) # should be sum over Nspin
    Poisson_solve_sum!(Ham.grid, ρ, Vhartree)
    Vxc[:,ispin] = calc_Vxc_1d(Ham.xc_calc, rhoe_old)
    Vtot[:,ispin] = Vhartree[:] + Vxc[:,ispin]

    A_old = calc_free_energy(Ham, Focc_matrix, hmat, Vhartree, rhoe_old)
    println("A_old = ", A_old)
    
    #for iter_inner in 1:2
    #    change_f, _ = inner_loop!(Ham, Focc_matrix, psi, hmat, NmaxIter_inner = 1)
    #    println("iter_inner = $(iter_inner) change_f = $(change_f)")
    #end

    change_f, _ = inner_loop!(Ham, Focc_matrix, psi, hmat, NmaxIter_inner = 2)
    println("change_f = $(change_f)")

    rhoe_new = similar(rhoe)
    edft_calc_rhoe!(Ham, Focc_matrix, psi, rhoe_new)
    #
    ρ = reshape(rhoe_new, Npoints) # should be sum over Nspin
    Poisson_solve_sum!(Ham.grid, ρ, Vhartree)
    Vxc[:,ispin] = calc_Vxc_1d(Ham.xc_calc, rhoe_new)
    
    A_new = calc_free_energy(Ham, Focc_matrix, hmat, Vhartree, rhoe_new)
    println("A_new = ", A_new)

    dn = density_difference(rhoe_new, rhoe_old, hx)
    println("dn = ", dn)

    @infiltrate

end

