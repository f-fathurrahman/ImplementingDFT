using Printf
using LinearAlgebra
using SparseArrays
using BenchmarkTools

include("FD3dGrid.jl")
include("build_nabla2_matrix.jl")
include("../common/supporting_functions.jl")

function test_multiplication()

    println()
    println("Benchmarking sparse matrix multiplication")
    println()

    Nx = 50
    Ny = 50
    Nz = 50
    grid = FD3dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny, (-5.0,5.0), Nz )

    println("Nx = ", Nx)
    println("Ny = ", Ny)
    println("Nz = ", Nz)

    print("Building ∇2 (1st call): ")
    @time ∇2 = build_nabla2_matrix( grid )
    
    print("Building ∇2 (2nd call): ")
    @time ∇2 = build_nabla2_matrix( grid )

    psi1 = rand( grid.Nx, grid.Ny, grid.Nz )
    psi2 = rand( grid.Npoints )

    println("3d array, using views:")
    @btime begin
       @views res1 = $∇2*$psi1[:]
    end

    println("3d array, not using views:")
    @btime begin
       res1 = $∇2*$psi1[:]
    end

    println("Using one-dimensional array")
    @btime begin
       res2 = $∇2*$psi2
    end

    println("Pass here")
end
test_multiplication()



