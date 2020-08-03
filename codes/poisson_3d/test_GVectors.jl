using Printf
using LinearAlgebra
using SparseArrays
using IncompleteLU
using AlgebraicMultigrid
using FFTW

include("INC_poisson_3d.jl")
include("../common/GVectors.jl")
include("Poisson_solve_fft.jl")

function test_main( NN::Array{Int64} )
    AA = [-8.0, -8.0, -8.0]
    BB = [ 8.0,  8.0,  8.0]
    grid = FD3dGrid( NN, AA, BB, pbc=(true,true,true) )
    println(grid)
    gvec = GVectors(grid)
    println(gvec)
end

test_main([14,14,14])