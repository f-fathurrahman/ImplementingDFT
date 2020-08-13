function KS_solve_SCF_potmix!(
    Ham::Hamiltonian, psi::Array{Float64,2};
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

    Rhoe = zeros(Float64,Npoints)
    Rhoe_new = zeros(Float64,Npoints)
    
    if guess_density == :random
        calc_rhoe!( Ham, psi, Rhoe )
        update!( Ham, Rhoe )
        evals = zeros(Float64,Nstates)
    else
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

    Rhoe_old = copy(Rhoe)
    Vxc_inp = zeros(Float64,Npoints)
    VHa_inp = zeros(Float64,Npoints)

    ethr_evals_last = 1e-5
    ethr = 0.1

    Ham.energies.NN = calc_E_NN( Ham.atoms, Ham.pspots )

    Nconverges = 0

    for iterSCF in 1:NiterMax

        # determine convergence criteria for diagonalization
        if iterSCF == 1
            ethr = 0.1
        elseif iterSCF == 2
            ethr = 0.01
        else
            ethr = ethr/5.0
            ethr = max( ethr, ethr_evals_last )
        end

        evals = diag_func( Ham, psi, Ham.precKin, tol=ethr,
                           Nstates_conv=Ham.electrons.Nstates_occ )
        if diag_func == diag_davidson!
            psi = psi*sqrt(dVol) # for diag_davidson
        else
            psi = psi/sqrt(dVol) # renormalize
        end

        if use_smearing
            E_f, Ham.energies.mTS =
            update_Focc!( Ham.electrons.Focc, smear_func, smear_func_entropy,
                      evals, Float64(Ham.electrons.Nelectrons), kT )
            @printf("Fermi energy = %18.10f\n", E_f)
        end

        calc_rhoe!( Ham, psi, Rhoe )

        # Save the old (input) potential, before updating the potential
        Vxc_inp[:] = Ham.V_XC
        VHa_inp[:] = Ham.V_Hartree

        update!( Ham, Rhoe )

        # Now Ham.potentials contains new (output) potential

        calc_energies!( Ham, psi )
        Etot = sum( Ham.energies )

        dRhoe = sum(abs.(Rhoe - Rhoe_old))/Npoints # MAE
        dEtot = abs(Etot - Etot_old)

        @printf("SCF_potmix: %5d %18.10f %18.10e %18.10e\n", iterSCF, Etot, dEtot, dRhoe)
        if Etot_old - Etot < 0.0
            println("WARNING: Energy is not reducing")
        end

        if dEtot < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        if Nconverges >= 2
            @printf("\nSCF_potmix is converged in iter: %d\n", iterSCF)
            @printf("\nOccupations and eigenvalues:\n")
            for i in 1:Nstates
                @printf("%3d %8.5f %18.10f\n", i, Ham.electrons.Focc[i], evals[i])
            end
            break
        end

        # Mix potential
        Ham.V_Hartree = betamix*Ham.V_Hartree + (1-betamix)*VHa_inp
        Ham.V_XC = betamix*Ham.V_XC + (1-betamix)*Vxc_inp

        Etot_old = Etot
        Rhoe_old = copy(Rhoe)
        flush(stdout)

    end

    println(Ham.energies)

    return
end