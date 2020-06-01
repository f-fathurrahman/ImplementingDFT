using Printf
using LinearAlgebra
using SparseArrays
using BenchmarkTools

include("../FD2d/FD2dGrid.jl")
include("../FD2d/build_nabla2_matrix.jl")
include("../common/supporting_functions.jl")

include("ILU0Preconditioner.jl")

function main()
    Nx = 5
    Ny = 4
    grid = FD2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )
    ∇2 = build_nabla2_matrix( grid, stencil_order=3 )

    prec = ILU0Preconditioner(-∇2)

    x = randn(Nx*Ny)
    println("dot(x,x) before = ", dot(x,x))

    ldiv!(prec, x) 

    println("dot(x,x) after  = ", dot(x,x))
end

main()