push!(LOAD_PATH, "../")

using Printf
using LinearAlgebra
using KSDFT1d
using Serialization

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
    E_NN = 2*2/2 # Z_i*Z_j/r_ij
    β_mix = 0.2
    
    Ekin = 0.0
    Ehartree = 0.0
    Eion = 0.0
    Exc = 0.0
    Etot = 0.0

    psi = zeros(Float64, Npoints, Nstates)

    for iter_scf in 1:100

        # Hamiltonian
        @views Vtot[:,1] = Vion[:] + Vhartree[:] + Vxc[:,1]
        @views Hmat[:,:] = Kmat + diagm( 0 => Vtot[:,1] )
    
        # Solve the eigenproblem
        evals_all, evecs_all = eigen( Hmat )   
        psi[:,:] .= evecs_all[:,1:Nstates]
        evals = evals_all[1:Nstates]

        # Renormalize
        psi[:] = psi[:]/sqrt(hx)
        
        @printf("Eigenvalues\n")
        for ist in 1:Nstates
            @printf("%5d %18.10f\n", ist, evals[ist])
        end

        calc_rhoe!(Ham, psi, rhoe_new)
        println("integ rhoe_new = ", sum(rhoe_new)*hx)

        epsxc[:] = calc_epsxc_1d(Ham.xc_calc, rhoe_new[:,1])
        Ekin = calc_E_kin(Ham, psi)
        Ehartree = 0.5*dot(rhoe_new[:,1], Vhartree)*hx
        Eion = dot(rhoe_new, Vion)*hx
        Exc = dot(rhoe_new, epsxc)*hx
        Etot = Ekin + Ehartree + Eion + Exc + E_NN

        ΔE = abs(Etot - Etot_old)
        mae_rhoe = sum(abs.(rhoe - rhoe_new))/Npoints
        @printf("%3d %18.10f %10.5e %10.5e\n", iter_scf, Etot, ΔE, mae_rhoe)

        if mae_rhoe < 1e-7
            println("Converged")
            break
        end

        # Mix
        if iter_scf >= 2
            rhoe[:] = β_mix*rhoe_new[:] + (1 - β_mix)*rhoe[:]
        else
            rhoe[:] = rhoe_new[:]
        end
        Etot_old = Etot

        # Update the potentials
        ρ = reshape(rhoe, Npoints)
        Poisson_solve_sum!(Ham.grid, ρ, Vhartree)
        Vxc[:] = calc_Vxc_1d(Ham.xc_calc, rhoe)

    end

    serialize("TEMP_psi.dat", psi)

    @printf("-----------------------------\n")
    @printf("Total energy components\n")
    @printf("-----------------------------\n")
    @printf("Ekin     = %18.10f\n", Ekin)
    @printf("Eion     = %18.10f\n", Eion)
    @printf("Ehartree = %18.10f\n", Ehartree)
    @printf("Exc      = %18.10f\n", Exc)
    @printf("E_NN     = %18.10f\n", E_NN)
    @printf("-----------------------------\n")
    @printf("Etot     = %18.10f\n", Etot)


end

main()
