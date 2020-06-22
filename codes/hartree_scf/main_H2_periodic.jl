push!(LOAD_PATH, pwd())

using Printf
using LinearAlgebra
using SparseArrays
using AlgebraicMultigrid
using Random

using MyModule

using FFTW
using SpecialFunctions: erf
include("../common/potential_H_atom.jl")

function main()

    Random.seed!(1234)

    AA = [ 0.0,  0.0,  0.0]
    BB = [16.0, 16.0, 16.0]
    NN = [41, 41, 41]

    grid = FD3dGrid( NN, AA, BB, pbc=(true,true,true) )

    atoms = Atoms( xyz_string=
        """
        2

        H   0.75  0.0  0.0
        H  -0.75  0.0  0.0
        """, in_bohr=true, pbc=(true,true,true))

    gvec = GVectors(grid)
    #V_Ps_loc = pot_H_atom_G(atoms, grid, gvec)
    V_Ps_loc = pot_Hps_GTH_G(atoms, grid, gvec)

    Nstates = 1
    Nelectrons = 2

    Ham = Hamiltonian(atoms, grid, V_Ps_loc, Nelectrons=Nelectrons)

    Npoints = grid.Npoints
    dVol = grid.dVol

    psi = rand(Float64,Npoints,Nstates)
    ortho_sqrt!(psi)
    psi = psi/sqrt(dVol)

    println("Test normalization")
    for i in 1:Nstates
        @printf("%18.10f\n", dot(psi[:,i], psi[:,i])*dVol )
    end

    Rhoe = calc_rhoe( Ham, psi )
    @printf("Integrated Rhoe = %18.10f\n", sum(Rhoe)*dVol)

    update!( Ham, Rhoe )

    evals = zeros(Float64,Nstates)
    Etot_old = 0.0
    dEtot = 0.0
    betamix = 0.5
    dRhoe = 0.0
    NiterMax = 100

    Ham.energies.NN = calc_E_NN( atoms.LatVecs, atoms, [1.0] )

    for iterSCF in 1:NiterMax

        evals = diag_LOBPCG!( Ham, psi, Ham.precKin, verbose_last=false )
        psi = psi/sqrt(dVol)

        Rhoe_new = calc_rhoe( Ham, psi )
        Rhoe = betamix*Rhoe_new + (1-betamix)*Rhoe

        update!( Ham, Rhoe )

        Etot = sum( calc_energies( Ham, psi ) )

        dRhoe = sum(abs.(Rhoe - Rhoe_new))/Npoints # MAE
        dEtot = abs(Etot - Etot_old)

        @printf("%5d %18.10f %18.10e %18.10e\n", iterSCF, Etot, dEtot, dRhoe)

        if dEtot < 1e-6
            @printf("Convergence is achieved in %d iterations\n", iterSCF)
            for i in 1:Nstates
                @printf("%3d %18.10f\n", i, evals[i])
            end
            break
        end

        Etot_old = Etot
    end

    println(Ham.energies)

end

main()
