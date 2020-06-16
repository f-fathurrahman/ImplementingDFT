push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

include("../common/potential_H_atom.jl")

function main()
    Random.seed!(1234)

    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [50, 50, 50]

    grid = FD3dGrid( NN, AA, BB )
    #grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    println(grid)

    atoms = Atoms(xyz_string="""
    1

    H   0.0   0.0   0.0
    """)

    V_Ps_loc = pot_Hps_HGH( atoms, grid )

    Nstates = 1
    Nelectrons = 1
    Ham = Hamiltonian( atoms, grid, V_Ps_loc, Nelectrons=1 )

    Npoints = grid.Npoints
    dVol = grid.dVol
    psi = rand(Float64,Npoints,Nstates)
    ortho_sqrt!(psi)
    psi = psi/sqrt(dVol)

    for i in 1:Nstates
        @printf("%18.10f\n", dot(psi[:,i], psi[:,i])*dVol )
    end

    Rhoe_new = zeros(Float64,Npoints)
    Rhoe = zeros(Float64,Npoints)

    calc_rhoe!( Ham, psi, Rhoe )
    @printf("Integrated Rhoe = %18.10f\n", sum(Rhoe)*dVol)

    update!( Ham, Rhoe )

    evals = zeros(Float64,Nstates)
    Etot_old = 0.0
    dEtot = 0.0
    betamix = 0.5
    dRhoe = 0.0
    NiterMax = 100

    for iterSCF in 1:NiterMax

        evals = diag_LOBPCG!( Ham, psi, Ham.precKin, verbose_last=false )
        #evals = diag_Emin_PCG!( Ham, psi, Ham.precKin, verbose_last=false )

        #psi = psi*sqrt(dVol) # for diag_davidson
        #evals = diag_davidson!( Ham, psi, Ham.precKin, verbose_last=false )

        psi = psi/sqrt(dVol) # renormalize

        #Rhoe_new = calc_rhoe( Ham, psi )
        calc_rhoe!( Ham, psi, Rhoe_new )

        Rhoe = betamix*Rhoe_new + (1-betamix)*Rhoe

        update!( Ham, Rhoe )

        calc_energies!( Ham, psi )
        Etot = sum( Ham.energies )

        dRhoe = sum(abs.(Rhoe - Rhoe_new))/Npoints
        dEtot = abs(Etot - Etot_old)

        @printf("%5d %18.10f %18.10e %18.10e\n", iterSCF, Etot, dEtot, dRhoe)

        if dEtot < 1e-6
            @printf("Convergence is achieved in %d iterations\n", iterSCF)
            @printf("\nEigenvalues:\n")
            for i in 1:Nstates
                @printf("%3d %18.10f\n", i, evals[i])
            end
            break
        end

        Etot_old = Etot
    end

    println(Ham.energies)

end

@time main()
