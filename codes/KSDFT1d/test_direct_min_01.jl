push!(LOAD_PATH, "./")

using Printf
using LinearAlgebra
using KSDFT1d

function create_atoms()
    Natoms = 2
    σ = ones(Float64, Natoms)*(1.0)
    masses = ones(Float64, Natoms)*42000.0
    Zvals = ones(Float64, Natoms)*2
    L = 10.0
    atpos = zeros(Float64, Natoms)
    atpos[1] = -1.0
    atpos[2] =  1.0
    return Atoms1d( atpos, Zvals, σ, masses, L )
end

function pot_gaussian( x, x0 )
    return -exp(-(x - x0)^2)/sqrt(1.0/π)^3
end

function init_Vions!(Ham)
    atpos = Ham.atoms.positions
    Natoms = Ham.atoms.Natoms
    for ia in 1:Natoms
        Ham.potentials.Ions[:] += pot_gaussian.(Ham.grid.x, atpos[ia])
    end
    return
end

function init_Hamiltonian()
    atoms = create_atoms()
    Ham = Hamiltonian1d(atoms, 51)
    init_Vions!(Ham)
    return Ham
end

# Evaluate gradient at psi
# Ham matrix is not updated
function calc_grad( Ham, psi )
    
    ispin = 1
    Focc = Ham.electrons.Focc
    hx = Ham.grid.hx

    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    #
    grad = zeros(Float64, Nbasis, Nstates)

    Hpsi = Ham*psi
    for i in 1:Nstates
        @views grad[:,i] .= Hpsi[:,i]
        for j in 1:Nstates
            @views λ = dot(psi[:,j], Hpsi[:,i]) * hx
            @views grad[:,i] .= grad[:,i] .-  λ * psi[:,j]
        end
        grad[:,i] .= Focc[i,ispin]*grad[:,i]
    end

    # the usual case of constant occupation numbers
    if all(Focc[:,ispin] .== 2.0) == true || all(Focc[:,ispin] .== 1.0) == true
        # immediate return
        return grad
    end

    #F = Matrix(Diagonal(Focc[:,ispin]))
    #ℍ = psi' * Hpsi * hx
    #HFH = ℍ*F - F*ℍ
    #ℚ = 0.5*HFH
    #@views grad[:,:] .+= psi*ℚ # additional contributions
    
    return grad
end

function ortho_sqrt( psi::Array{Float64,2} )
    Udagger = inv(sqrt(psi'*psi))
    return psi*Udagger
end

function ortho_sqrt!( psi::Array{Float64,2} )
    Udagger = inv(sqrt(psi'*psi))
    psi[:,:] = psi*Udagger
    return
end


# Ham will be modified (the potential terms)
function calc_KohnSham_Etotal!(Ham, psi)
    hx = Ham.grid.hx    
    Npoints = Ham.grid.Npoints
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    Vtot = Ham.potentials.Total
    rhoe = Ham.rhoe

    # Electron density
    calc_rhoe!(Ham, psi, rhoe)
    println("integ rhoe = ", sum(rhoe)*hx)

    # Update the potentials
    ρ = reshape(rhoe, Npoints)
    Poisson_solve_sum!(Ham.grid, ρ, Vhartree)
    Vxc[:] = calc_Vxc_1d(Ham.xc_calc, rhoe)
    @views Vtot[:,1] = Vion[:] + Vhartree[:] + Vxc[:,1]

    epsxc = calc_epsxc_1d(Ham.xc_calc, rhoe[:,1])
    Ekin = calc_E_kin(Ham, psi)
    Ehartree = 0.5*dot(rhoe[:,1], Vhartree)*hx
    Eion = dot(rhoe, Vion)*hx
    Exc = dot(rhoe, epsxc)*hx
    Eelec = Ekin + Ehartree + Eion + Exc
    # only electronic parts, E_NN should be added separately

    @printf("-----------------------------\n")
    @printf("Total energy components\n")
    @printf("-----------------------------\n")
    @printf("Ekin     = %18.10f\n", Ekin)
    @printf("Eion     = %18.10f\n", Eion)
    @printf("Ehartree = %18.10f\n", Ehartree)
    @printf("Exc      = %18.10f\n", Exc)
    @printf("-----------------------------\n")
    @printf("Eelec     = %18.10f\n", Eelec)

    return Eelec
end

function prec_invK(Ham::Hamiltonian1d, v)
    return inv(Ham.Kmat)*v
end

function main()

    Ham = init_Hamiltonian()

    hx = Ham.grid.hx
    Npoints = Ham.grid.Npoints
    Nelectrons = Ham.electrons.Nelectrons
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc

    Kmat = Ham.Kmat
    Vtot = Ham.potentials.Total
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    rhoe = Ham.rhoe

    Nspin = 1
    Hmat = zeros(Float64, Npoints, Npoints)
    epsxc = zeros(Float64, Npoints)
    rhoe_new = zeros(Float64, Npoints, Nspin)

    Etot = Inf
    Etot_old = Etot
    E_NN = calc_E_NN(Ham.atoms)

    psi = ortho_sqrt( rand(Float64, Npoints, Nstates) )
    psi[:] = psi[:]/sqrt(hx)
    display(psi' * psi * hx)


    α_t = 3e-5
    g = zeros(Float64, Npoints, Nstates)
    Kg = zeros(Float64, Npoints, Nstates)
    gt = zeros(Float64, Npoints, Nstates)
    d = zeros(Float64, Npoints, Nstates)
    psic = zeros(Float64, Npoints, Nstates)

    for iterEmin in 1:50
        
        Eelec = calc_KohnSham_Etotal!(Ham, psi)
        Etot = Eelec + E_NN

        println("iterEmin = ", iterEmin, " Etot = ", Etot)

        g[:,:] .= calc_grad(Ham, psi)
        Kg[:,:] .= prec_invK(Ham, g)

        println("dot(g,g) = ", dot(g,g)*hx)

        @views d[:] .= -Kg[:]
        psic[:] = ortho_sqrt( psi + α_t*d )
        psic[:] = psic[:]/sqrt(hx)

        _ = calc_KohnSham_Etotal!(Ham, psic)
        gt[:,:] = calc_grad(Ham, psic)

        denum = real(sum(conj(g-gt).*d))
        if denum != 0.0
            α = abs( α_t*real(sum(conj(g).*d))/denum )
        else
            α = 0.0
        end
        #denum = sum( (g - gt) .* d )
        #if denum != 0.0
        #    α = abs( α_t * sum(g .* d)/denum )
        #else
        #    α = 0.0
        #end
        println("α = ", α)

        # Update wavefunction
        psi[:,:] .= ortho_sqrt(psi + α*d)
        psi[:] = psi[:]/sqrt(hx)
    end
    println("E_NN = ", E_NN)

    Hr = psi' * (Ham*psi) * hx
    band_energies = eigvals(Hr)
    display(band_energies); println()


end

main()
