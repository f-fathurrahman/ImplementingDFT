using Printf
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using IncompleteLU
using AlgebraicMultigrid
using Random

include("../../3d/FD3dGrid.jl")
include("../../3d/build_nabla2_matrix.jl")
include("../../diag_Emin_PCG.jl")
include("../../diag_davidson.jl")
include("../../diag_LOBPCG.jl")
include("../../ortho_sqrt.jl")
include("../../supporting_functions.jl")
include("../../3d_poisson/Poisson_solve_PCG.jl")
include("calc_rhoe.jl")
include("Hamiltonian.jl")
include("calc_energies.jl")

function pot_harmonic( fdgrid::FD3dGrid; ω=1.0, center=[0.0, 0.0, 0.0] )
    Npoints = fdgrid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = fdgrid.r[1,i] - center[1]
        y = fdgrid.r[2,i] - center[2]
        z = fdgrid.r[3,i] - center[3]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 + z^2 )
    end
    return Vpot
end

function main()

    Random.seed!(1234)
    
    AA = [0.0, 0.0, 0.0]
    BB = [6.0, 6.0, 6.0]
    NN = [25, 25, 25]
    
    fdgrid = FD3dGrid( NN, AA, BB )

    # follow the potential used in Arias DFT++ tutorial
    my_pot_harmonic( fdgrid ) = pot_harmonic( fdgrid, ω=2, center=[3.0, 3.0, 3.0] )

    Ham = Hamiltonian( fdgrid, my_pot_harmonic )

    Nbasis = prod(NN)

    dVol = fdgrid.dVol

    Nstates = 4
    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi)
    psi = psi/sqrt(dVol)

    for i in 1:Nstates
        @printf("%18.10f\n", dot(psi[:,i], psi[:,i])*dVol )
    end

    Rhoe = calc_rhoe( psi )
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
        psi = psi/sqrt(dVol)

        Rhoe_new = calc_rhoe( psi )

        Rhoe = betamix*Rhoe_new + (1-betamix)*Rhoe

        update!( Ham, Rhoe )
        
        Etot = sum( calc_energies( Ham, psi ) )
        
        dRhoe = norm(Rhoe - Rhoe_new)
        dEtot = abs(Etot - Etot_old)

        @printf("%5d %18.10f %18.10e %18.10e\n", iterSCF, Etot, dEtot, dRhoe)

        if dEtot < 1e-6
            @printf("Convergence is achieved in %d iterations\n", iterSCF)
            break
        end

        Etot_old = Etot
    end
end

main()



