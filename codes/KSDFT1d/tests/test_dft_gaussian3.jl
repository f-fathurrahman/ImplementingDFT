push!(LOAD_PATH, "../")

using Printf
using LinearAlgebra
using Serialization

using KSDFT1d

include("../BroydenMixer.jl")

function create_atoms()
    Natoms = 3
    σ = ones(Float64, Natoms)*(1.0)
    masses = ones(Float64, Natoms)*42000.0
    Zvals = ones(Float64, Natoms)*5.0
    L = 10.0
    atpos = zeros(Float64, Natoms)
    atpos[1] = -1.5
    atpos[2] =  0.1
    atpos[3] =  1.5
    return Atoms1d( atpos, Zvals, σ, masses, L )
end

function pot_gaussian( x, x0 )
    return -25.0*exp(-4.5*(x-x0)^2)
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
    Ham = Hamiltonian1d(atoms, 51, Nstates_extra=6)
    init_Vions!(Ham)
    return Ham
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

    psi = zeros(Float64, Npoints, Nstates)

    Etot = Inf
    Etot_old = Etot
    Ham.energies.NN = calc_E_NN(Ham.atoms)

    E_NN = Ham.energies.NN

    betamix = 0.5
    mixer = BroydenMixer(rhoe, betamix, mixdim=8)

    Focc = Ham.electrons.Focc
    use_smearing = true
    kT = 0.1*eV2Ha #0.1 # 0.3*eV2Ha

    println("kT = ", kT)

    evals = Ham.electrons.ebands

    for iter_scf in 1:200

        # Hamiltonian
        @views Vtot[:,1] = Vion[:] + Vhartree[:] + Vxc[:,1]
        @views Hmat[:,:] = Kmat + diagm( 0 => Vtot[:,1] )
    
        # Solve the eigenproblem
        evals_all, evecs_all = eigen( Hmat )   
        psi[:,:] .= evecs_all[:,1:Nstates]
        evals[:,1] .= evals_all[1:Nstates]

        # Renormalize
        psi[:] = psi[:]/sqrt(hx)

        if use_smearing
            E_f, mTS =
            update_Focc!( Focc, smear_fermi, smear_fermi_entropy,
                      evals, Float64(Nelectrons), kT )
            @printf("Fermi energy = %18.10f\n", E_f)
        end

        @printf("Eigenvalues\n")
        for ist in 1:Nstates
            @printf("%5d %18.10f occ=%10.5f\n", ist, evals[ist,1], Focc[ist,1])
        end

        calc_rhoe!(Ham, psi, rhoe_new)
        println("integ rhoe_new = ", sum(rhoe_new)*hx)

        Ekin = calc_E_kin(Ham, psi)
        Ehartree = 0.5*dot(rhoe_new[:,1], Vhartree)*hx
        Eion = dot(rhoe_new, Vion)*hx

        epsxc[:] = calc_epsxc_1d(Ham.xc_calc, rhoe_new[:,1])
        Exc = dot(rhoe_new, epsxc)*hx
        
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
        println("integ rhoe after mix = ", sum(rhoe)*hx)

        Etot_old = Etot

        # Update the potentials
        ρ = reshape(rhoe, Npoints)
        Poisson_solve_sum!(Ham.grid, ρ, Vhartree)
        Vxc[:] = calc_Vxc_1d(Ham.xc_calc, rhoe)

    end

    serialize("TEMP_psi.dat", psi)
    serialize("TEMP_evals.dat", evals)

    println(Ham.energies)

end

main()
