# XXX rhoe is taken from Ham.rhoe so it might converges early if called
# multiple times consecutively
#
function solve_scf!(Ham)

    dx = Ham.grid.dx
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

    psi = zeros(Float64, Npoints, Nstates)

    Etot = Inf
    Etot_old = Etot
    Ham.energies.NN = calc_E_NN(Ham.atoms)

    E_NN = Ham.energies.NN

    betamix = 0.5
    mixer = BroydenMixer(rhoe, betamix, mixdim=8)

    Focc = Ham.electrons.Focc
    use_smearing = Ham.electrons.use_smearing
    kT = Ham.electrons.kT

    evals = Ham.electrons.ebands

    Nconverges = 0

    for iter_scf in 1:200

        println("\nBegin iter_scf = ", iter_scf)

        # Hamiltonian
        @views Vtot[:,1] = Vion[:] + Vhartree[:] + Vxc[:,1]
        @views Hmat[:,:] = Kmat + diagm( 0 => Vtot[:,1] )
    
        # Solve the eigenproblem
        evals_all, evecs_all = eigen( Hmat )   
        psi[:,:] .= evecs_all[:,1:Nstates]
        evals[:,1] .= evals_all[1:Nstates]

        # Renormalize
        psi[:] = psi[:]/sqrt(dx)

        if use_smearing
            Ham.electrons.E_fermi, mTS =
            update_Focc!( Focc, smear_fermi, smear_fermi_entropy,
                      evals, Float64(Nelectrons), kT )
        end

        calc_rhoe!(Ham, psi, rhoe_new)
        println("integ rhoe_new = ", sum(rhoe_new)*dx)

        Ekin = calc_E_kin(Ham, psi)
        Ehartree = 0.5*dot(rhoe_new[:,1], Vhartree)*dx
        Eion = dot(rhoe_new, Vion)*dx

        epsxc[:] = calc_epsxc_1d(Ham.xc_calc, rhoe_new[:,1])
        Exc = dot(rhoe_new, epsxc)*dx
        
        # Set the internal variables
        Ham.energies.Kinetic = Ekin
        Ham.energies.Hartree = Ehartree
        Ham.energies.Ion = Eion
        Ham.energies.XC = Exc
        Ham.energies.mTS = mTS

        Etot = Ekin + Ehartree + Eion + Exc + mTS + E_NN

        ΔE = abs(Etot - Etot_old)
        mae_rhoe = sum(abs.(rhoe - rhoe_new))/Npoints
        @printf("%3d %18.10f %10.5e %10.5e\n", iter_scf, Etot, ΔE, mae_rhoe)

        if mae_rhoe < 1e-7
            Nconverges += 1
        else
            Nconverges = 0
        end

        if Nconverges >= 2
            println("Converged: mae_rhoe = ", mae_rhoe, " ΔE = ", ΔE)
            break
        end

        # Mix
        if iter_scf >= 2
            #rhoe[:] = betamix*rhoe_new[:] + (1 - betamix)*rhoe[:]
            do_mix!(mixer, rhoe, rhoe_new, iter_scf)
        else
            rhoe[:] = rhoe_new[:]
        end
        println("integ rhoe after mix = ", sum(rhoe)*dx)

        Etot_old = Etot

        # Update the potentials
        ρ = reshape(rhoe, Npoints)
        Poisson_solve_sum!(Ham.grid, ρ, Vhartree)
        Vxc[:] = calc_Vxc_1d(Ham.xc_calc, rhoe)

    end

    @printf("Eigenvalues\n")
    for ist in 1:Nstates
        @printf("%5d %18.10f occ=%10.5f\n", ist, evals[ist,1], Focc[ist,1])
    end
    @printf("Fermi energy = %18.10f\n", Ham.electrons.E_fermi)

    serialize("TEMP_scf_psi.dat", psi)
    serialize("TEMP_scf_evals.dat", evals)

    println(Ham.energies)

end


#main()
