function calc_grad!( Ham, ψ, g, Hsub )
    Nstates = size(ψ,2)
    Focc = Ham.electrons.Focc
    Hψ = op_H( Ham, ψ )
    Hsub[:] = ψ' * Hψ * Ham.grid.dVol
    Hψ = Hψ - ψ*Hsub
    for ist in 1:Nstates
        g[:,ist] = Focc[ist] * Hψ[:,ist]
    end
    return
end

function calc_grad!( Ham, ψ, g )
    Nstates = size(ψ,2)
    Focc = Ham.electrons.Focc
    Hψ = op_H( Ham, ψ )
    Hsub = ψ' * Hψ * Ham.grid.dVol
    Hψ = Hψ - ψ*Hsub
    for ist in 1:Nstates
        g[:,ist] = Focc[ist] * Hψ[:,ist]
    end
    #println(dot(g,g))
    return
end

function calc_energies_grad!( Ham, psi, Rhoe, g, Kg, Hsub )
    calc_rhoe!( Ham, psi, Rhoe )
    update!( Ham, Rhoe )
    calc_energies!( Ham, psi )
    calc_grad!( Ham, psi, g, Hsub )
    # apply preconditioner
    Kg[:,:] = g[:,:]
    Nstates = size(psi,2)
    for ist in 1:Nstates
        @views ldiv!(Ham.precKin, Kg[:,ist])
    end
    return sum( Ham.energies )
end


function calc_energies_only!( Ham, psi, Rhoe )
    calc_rhoe!( Ham, psi, Rhoe )
    update!( Ham, Rhoe )
    calc_energies!( Ham, psi )
    return sum( Ham.energies )
end


function linmin_grad!( Ham, psi, g, d; αt = 3e-5 )

    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    
    psic = zeros(Float64,Nbasis,Nstates)
    gt = zeros(Float64,Nbasis,Nstates)
    
    dVol= Ham.grid.dVol

    psic = psi + αt*d
    ortho_sqrt!(psic)
    psic = psic/sqrt(dVol)

    Rhoe = calc_rhoe( Ham, psic )
    update!( Ham, Rhoe )
    calc_grad!( Ham, psic, gt )

    denum = dot(g - gt, d)*dVol

    if denum != 0.0
        α = abs( αt * dot(g, d)*dVol/denum )
    else
        α = 0.0
    end
    return α
end

function constrain_search_dir!( d, psi, dVol )
    d[:,:] = d[:,:] - psi * ( psi' * d * dVol )
    return
end

function KS_solve_Emin_PCG!(
    Ham::Hamiltonian, psi::Array{Float64,2};
    α_t=3e-5, NiterMax=200,
    i_cg_beta=2, etot_conv_thr=1e-6
)
    Nbasis = size(psi,1)
    Nstates = size(psi,2)

    dVol = Ham.grid.dVol

    g = zeros(Float64,Nbasis,Nstates)
    Kg = zeros(Float64,Nbasis,Nstates)
    gPrev = zeros(Float64,Nbasis,Nstates)
    d = zeros(Float64,Nbasis,Nstates)

    Rhoe = zeros(Float64,Nbasis)
    Hsub = zeros(Nstates,Nstates)
    Rhoe_old = zeros(Float64,Nbasis)

    Ham.energies.NN = calc_E_NN( Ham.atoms, Ham.pspots )

    Etot = calc_energies_grad!( Ham, psi, Rhoe, g, Kg, Hsub )

    d[:,:] = -Kg[:,:]

    constrain_search_dir!( d, psi, dVol )

    gPrevUsed = true

    α = 0.0
    β = 0.0
    gKnorm = 0.0
    gKnormPrev = 0.0
    force_grad_dir = true

    Etot_old = Etot
    Nconverges = 0
    
    cg_test = 0.0

    @printf("\n")
    @printf("Minimizing Kohn-Sham energy using PCG\n")
    @printf("-------------------------------------\n")
    @printf("NiterMax  = %d\n", NiterMax)
    @printf("α_t       = %e\n", α_t)
    @printf("conv_thr  = %e\n", etot_conv_thr)
    if i_cg_beta == 1
        @printf("Using Fletcher-Reeves formula for CG_BETA\n")
    elseif i_cg_beta == 2
        @printf("Using Polak-Ribiere formula for CG_BETA\n")
    elseif i_cg_beta == 3
        @printf("Using Hestenes-Stiefeld formula for CG_BETA\n")
    else
        @printf("Using Dai-Yuan formula for CG_BETA\n")
    end
    @printf("\n")


    for iter in 1:NiterMax

        gKnorm = dot(g,Kg)*dVol
        if !force_grad_dir
            dotgd = dot(g,d)*dVol
            if gPrevUsed
                dotgPrevKg = dot(gPrev, Kg)*dVol
            else
                dotgPrevKg = 0.0
            end
            β = (gKnorm - dotgPrevKg)/gKnormPrev # Polak-Ribiere
        end

        if β < 0.0
            #println("Resetting β")
            β = 0.0
        end

        force_grad_dir = false
        
        if gPrevUsed
            gPrev[:,:] = g[:,:]
        end
        gKnormPrev = gKnorm

        # Update search direction
        d[:] = -Kg[:] + β*d[:]

        constrain_search_dir!( d, psi, dVol )

        α = linmin_grad!( Ham, psi, g, d )

        Rhoe_old = copy(Ham.rhoe)

        # Update psi
        psi[:] = psi[:] + α*d[:]
        #display(psi'*psi); println()
        #exit()
        #ortho_sqrt!(psi, dVol)
        ortho_sqrt!(psi)
        psi = psi/sqrt(dVol)

        calc_rhoe!(Ham, psi, Rhoe)
        update!(Ham, Rhoe)

        # Calc new Etot and gradient
        Etot = calc_energies_grad!( Ham, psi, Rhoe, g, Kg, Hsub )
        
        diffE = Etot_old - Etot
        @printf("Emin_PCG step %8d = %18.10f  %12.7e\n", iter, Etot, diffE)
        if diffE < 0.0
            println("*** WARNING: Etot is not decreasing")
        end

        if abs(diffE) < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        if Nconverges >= 2
            @printf("\nEmin_PCG_dot is converged in iter: %d\n", iter)
            break
        end

        Etot_old = Etot

    end

    # Calculate eigenvalues
    evals, evecs = eigen(Symmetric(Hsub))
    psi = psi*evecs

    @printf("\n")
    @printf("----------------------------\n")
    @printf("Final Kohn-Sham eigenvalues:\n")
    @printf("----------------------------\n")
    @printf("\n")
    for i in 1:Nstates
        @printf("%3d %18.10f\n", i, evals[i])
    end

    @printf("\n")
    @printf("----------------------------\n")
    @printf("Final Kohn-Sham energies:\n")
    @printf("----------------------------\n")
    @printf("\n")
    println(Ham.energies, banner=false)

    return

end
