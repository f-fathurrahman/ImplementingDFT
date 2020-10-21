function KS_solve_SCF!(
    Ham::Hamiltonian, psis::Vector{Array{Float64,2}};
    NiterMax=200, betamix=0.5,
    etot_conv_thr=1e-6,
    diag_func=diag_LOBPCG!,
    use_smearing=false,
    smear_func=smear_fermi, smear_func_entropy=smear_fermi_entropy,
    kT=0.01,
    guess_density=:random
)

    Npoints = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol
    Nspin = Ham.Nspin

    Rhoe = zeros(Float64,Npoints,Nspin)
    Rhoe_new = zeros(Float64,Npoints,Nspin)
    
    Rhoe_tot = zeros(Float64,Npoints)
    Rhoe_new_tot = zeros(Float64,Npoints)
    magn_new = zeros(Float64,Npoints)
    magn = zeros(Float64,Npoints)
    
    if guess_density == :random
        calc_rhoe!( Ham, psis, Rhoe )
        update!( Ham, Rhoe )
        evals = zeros(Float64, Nstates, Nspin)
        println("Integ Rhoe = ", sum(Rhoe)*dVol)
    else
        error("Not supported yet")
        gen_gaussian_density!( Ham.grid, Ham.atoms, Ham.pspots, Rhoe )
        update!( Ham, Rhoe )
        #
        evals = diag_func( Ham, psi, Ham.precKin, tol=1e-3,
                           Nstates_conv=Ham.electrons.Nstates_occ )
        if diag_func == diag_davidson!
            psi = psi*sqrt(dVol) # for diag_davidson
        else
            psi = psi/sqrt(dVol) # renormalize
        end
        #
        calc_rhoe!(Ham, psi, Rhoe)
        update!( Ham, Rhoe )
    end

    Etot_old = 0.0
    dEtot = 0.0
    dRhoe = 0.0
    NiterMax = 100

    ethr_evals_last=1e-5
    ethr = 0.1

    # Mixer
    #mixer = LinearMixer(betamix)
    mixer = RPulayMixer( Npoints, 4, betamix; Nspin=Nspin)
    #mixer = AdaptiveLinearMixer(Npoints, betamix, 0.8, Nspin=Nspin)

    Ham.energies.NN = calc_E_NN( Ham.atoms, Ham.pspots )

    Nconverges = 0

    SMALL = eps()

    for iterSCF in 1:NiterMax

        println("\nBegin SCF iter: ", iterSCF)
        println("----------------------------")

        # determine convergence criteria for diagonalization
        if iterSCF == 1
            ethr = 0.1
        elseif iterSCF == 2
            ethr = 0.01
        else
            ethr = ethr/5.0
            ethr = max( ethr, ethr_evals_last )
        end

        for ispin in 1:Nspin
            Ham.ispin = ispin
            evals[:,ispin] = diag_func( Ham, psis[ispin], Ham.precKin, tol=ethr,
                                        Nstates_conv=Ham.electrons.Nstates_occ)
            if diag_func == diag_davidson!
                psis[ispin] = psis[ispin]*sqrt(dVol) # for diag_davidson
            else
                psis[ispin] = psis[ispin]/sqrt(dVol) # renormalize
            end
        end

        if use_smearing
            E_f, Ham.energies.mTS =
            update_Focc!( Ham.electrons.Focc, smear_func, smear_func_entropy,
                      evals, Float64(Ham.electrons.Nelectrons), kT )
            @printf("Fermi energy = %18.10f\n", E_f)
        end
        @printf("Nelectrons = %d\n", Ham.electrons.Nelectrons)
        @printf("sum(Focc)  = %18.10f\n", sum(Ham.electrons.Focc))

        @printf("\nOccupations and eigenvalues:\n")
        for ist in 1:Nstates
            @printf("%3d | %8.5f %18.10f", ist, Ham.electrons.Focc[ist,1], evals[ist,1])
            if Nspin == 2
                @printf(" | %8.5f %18.10f\n", Ham.electrons.Focc[ist,2], evals[ist,2])
            else
                @printf("\n")
            end
        end

        calc_rhoe!( Ham, psis, Rhoe_new )
        
        #integ_rho = sum(Rhoe_new)*dVol
        #@printf("integ_rho (before renormalized) = %18.10f\n", integ_rho)
        #for ip in 1:length(Rhoe_new)
        #    Rhoe_new[ip] = Ham.electrons.Nelectrons/integ_rho * Rhoe_new[ip]
        #end
        #@printf("integ Rhoe before mix = %18.10f\n", sum(Rhoe_new)*dVol)


        #println("Linear mixing: betamix = ", betamix)
        #Rhoe = betamix*Rhoe_new + (1-betamix)*Rhoe
        
        #println("Adaptive mixing: betamix = ", betamix)
        #mix_adaptive!( Rhoe, Rhoe_new, betamix, betav, df )
        
        mixer.iter = iterSCF
        do_mix!(mixer, Rhoe, Rhoe_new)

        #if Nspin == 2
        #    Rhoe_tot[:] = dropdims(sum(Rhoe,dims=2),dims=2)
        #    Rhoe_new_tot[:] = dropdims(sum(Rhoe_new,dims=2),dims=2)
        #    #
        #    magn[:] = Rhoe[:,1] -Rhoe[:,2]
        #    magn_new[:] = Rhoe_new[:,1] - Rhoe_new[:,2]
        #    #
        #    # Mix
        #    Rhoe_tot = betamix*Rhoe_new_tot + (1-betamix)*Rhoe_tot
        #    magn = betamix*magn_new + (1-betamix)*magn
        #    #
        #    Rhoe[:,1] = 0.5*(Rhoe_tot + magn)
        #    Rhoe[:,2] = 0.5*(Rhoe_tot - magn)
        #end

        for ip in 1:Npoints*Nspin
            if Rhoe[ip] < SMALL
                #println("Negative Rhoe is detected. Setting to SMALL")
                Rhoe[ip] = SMALL
            end
        end

        #integ_rho = sum(Rhoe)*dVol
        #@printf("integ Rhoe after mix (before renormalized) = %18.10f\n", integ_rho)
        #println("diff old and after mix: ", sum(Rhoe - Rhoe_new))
        #for ip in 1:length(Rhoe)
        #    Rhoe[ip] = Ham.electrons.Nelectrons/integ_rho * Rhoe[ip]
        #end
        #integ_rho = sum(Rhoe)*dVol
        #@printf("integ Rhoe after mix (after renormalized)  = %18.10f\n", integ_rho)

        if Nspin == 2
            @views smagn = sum(Rhoe[:,1] .- Rhoe[:,2])*dVol
            println("Integ magn = ", smagn)
        end

        update!( Ham, Rhoe )

        calc_energies!( Ham, psis )
        Etot = sum( Ham.energies )

        dRhoe = sum(abs.(Rhoe - Rhoe_new))/Npoints # MAE
        dEtot = abs(Etot - Etot_old)

        @printf("SCF: %5d %18.10f %18.10e %18.10e\n", iterSCF, Etot, dEtot, dRhoe)

        if dEtot < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        if Nconverges >= 2
            @printf("\nSCF is converged in iter: %d\n", iterSCF)
            break
        end

        Etot_old = Etot
    end

    println(Ham.energies)

    return
end