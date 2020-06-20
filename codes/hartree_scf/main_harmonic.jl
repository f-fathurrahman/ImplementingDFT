push!(LOAD_PATH, pwd())

using Printf
using LinearAlgebra
using SparseArrays
using AlgebraicMultigrid
using Random

#include("INC_hartree_scf.jl")
using MyModule

function pot_harmonic( grid; ω=1.0, center=[0.0, 0.0, 0.0] )
    Npoints = grid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = grid.r[1,i] - center[1]
        y = grid.r[2,i] - center[2]
        z = grid.r[3,i] - center[3]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 + z^2 )
    end
    return Vpot
end

function main()

    Random.seed!(1234)

    AA = [-3.0, -3.0, -3.0]
    BB = [3.0, 3.0, 3.0]
    NN = [25, 25, 25]

    grid = FD3dGrid( NN, AA, BB )

    V_Ps_loc = pot_harmonic( grid, ω=2 )

    Nelectrons = 8
    Nstates = round(Int64,Nelectrons/2)

    Ham = Hamiltonian( Atoms(), grid, V_Ps_loc, Nelectrons=Nelectrons )

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

    for iterSCF in 1:NiterMax

        evals = diag_LOBPCG!( Ham, psi, Ham.precKin, verbose_last=false )
        #evals = diag_Emin_PCG!( Ham, psi, Ham.precKin, verbose_last=true )
        psi = psi/sqrt(dVol)

        Rhoe_new = calc_rhoe( Ham, psi )
        #@printf("Integrated Rhoe_new = %18.10f\n", sum(Rhoe_new)*dVol)

        Rhoe = betamix*Rhoe_new + (1-betamix)*Rhoe
        #@printf("Integrated Rhoe     = %18.10f\n", sum(Rhoe)*dVol)

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
end

main()
