using Printf
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using IncompleteLU
using AlgebraicMultigrid
using Random
using Serialization

include("INC_sch_3d.jl")

function pot_harmonic( grid::FD3dGrid; ω=1.0 )
    Npoints = grid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = grid.r[1,i]
        y = grid.r[2,i]
        z = grid.r[3,i]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 + z^2 )
    end
    return Vpot
end

function main()

    Random.seed!(1234)

    Nx = 50
    Ny = 50
    Nz = 50
    grid = FD3dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny, (-5.0,5.0), Nz )

    ∇2 = build_nabla2_matrix( grid, stencil_order=9 )

    Vpot = pot_harmonic( grid )
    
    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    # may choose between these two
    println("Building preconditioner")
    
    #prec = ilu(-0.5*∇2)
    #prec = ilu(Ham) # this should result in faster convergence
    
    prec = aspreconditioner(ruge_stuben(Ham))
    #prec = aspreconditioner(smoothed_aggregation(Ham))
    #prec = aspreconditioner(ruge_stuben(-0.5*∇2))
    #prec = aspreconditioner(smoothed_aggregation(-0.5*∇2))
    
    #prec = ILU0Preconditioner(Ham)

    println("Done building preconditioner")

    Nstates = 10
    Npoints = Nx*Ny*Nz
    X = rand(Float64, Npoints, Nstates)
    ortho_sqrt!(X)
    
    #evals = diag_Emin_PCG!( Ham, X, prec, verbose=true )
    #evals = diag_davidson!( Ham, X, prec, verbose=true )
    evals = diag_LOBPCG!( Ham, X, prec, verbose=true )

    @printf("\n\nEigenvalues\n")
    for i in 1:Nstates
        @printf("%5d %18.10f\n", i, evals[i])
    end

    #serialize("wavefunc.data", X)

end

main()



