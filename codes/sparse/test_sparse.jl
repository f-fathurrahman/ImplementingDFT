using Printf
using LinearAlgebra
using SparseArrays
using BenchmarkTools

include("../FD2d/FD2dGrid.jl")
include("../FD2d/build_nabla2_matrix.jl")
include("../common/supporting_functions.jl")

include("ilu0.jl")

function main()
    Nx = 5
    Ny = 4
    grid = FD2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )
    ∇2 = build_nabla2_matrix( grid, stencil_order=3 )

    n = ∇2.n
    a = -∇2.nzval
    ja = ∇2.rowval
    ia = ∇2.colptr

    Nnzval = length(a)
    alu_ilu0 = zeros( Float64, Nnzval+1 )
    jlu_ilu0 = zeros( Int64, Nnzval+1 )
    ju_ilu0 = zeros( Int64, n )
    iw_ilu0 = zeros( Int64, n )

    init_ilu0!( n, a, ja, ia, alu_ilu0, jlu_ilu0, ju_ilu0, iw_ilu0 )

    y = randn(Nnzval)
    println("dot y before: ", dot(y,y))
    
    x = zeros(Nnzval)

    lusol!(n, y, x, alu_ilu0, jlu_ilu0, ju_ilu0)    

    println("dot x x = ", dot(x,x))

    println("Pass here")
end

main()