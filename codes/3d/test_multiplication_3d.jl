using LinearAlgebra
using SparseArrays
using BenchmarkTools

include("FD3dGrid.jl")
include("build_nabla2_matrix.jl")
include("../supporting_functions.jl")

function test_multiplication()
    Nx = 50
    Ny = 50
    Nz = 50
    fdgrid = FD3dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny, (-5.0,5.0), Nz )

    println("Nx = ", Nx)
    println("Ny = ", Ny)
    println("Nz = ", Nz)

    @time ∇2 = build_nabla2_matrix( fdgrid )
    @time ∇2 = build_nabla2_matrix( fdgrid )

    psi1 = rand( fdgrid.Nx, fdgrid.Ny, fdgrid.Nz )
    psi2 = rand( fdgrid.Npoints )

    println("Using views:")
    @btime begin
       @views res1 = $∇2*$psi1[:]
    end

    println("Not using views:")
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



