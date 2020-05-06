using Printf
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using IncompleteLU
using AlgebraicMultigrid
using Random
using Serialization

include("FD3dGrid.jl")
include("build_nabla2_matrix.jl")
include("../diag_Emin_PCG.jl")
include("../diag_davidson.jl")
include("../diag_LOBPCG.jl")
include("../ortho_sqrt.jl")
include("../supporting_functions.jl")

function pot_harmonic( fdgrid::FD3dGrid; ω=1.0 )
    Npoints = fdgrid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = fdgrid.r[1,i]
        y = fdgrid.r[2,i]
        z = fdgrid.r[3,i]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 + z^2 )
    end
    return Vpot
end

function main()

    Random.seed!(1234)

    Nx = 50
    Ny = 50
    Nz = 50
    fdgrid = FD3dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny, (-5.0,5.0), Nz )

    ∇2 = build_nabla2_matrix( fdgrid, func_1d=build_D2_matrix_5pt )

    Vpot = pot_harmonic( fdgrid )
    
    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    # may choose between these two
    println("Building preconditioner")
    #prec = ilu(-0.5*∇2)
    #prec = ilu(Ham) # this should result in faster convergence
    prec = aspreconditioner(ruge_stuben(Ham))
    #prec = aspreconditioner(smoothed_aggregation(Ham))
    #prec = aspreconditioner(ruge_stuben(-0.5*∇2))
    #prec = aspreconditioner(smoothed_aggregation(-0.5*∇2))
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

    serialize("wavefunc.data", X)

end

main()



