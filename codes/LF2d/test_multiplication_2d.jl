using LinearAlgebra
using SparseArrays
using BenchmarkTools
using Printf

include("../LF1d/init_LF1d_c_grid.jl")
include("../common/supporting_functions.jl")
include("LF2dGrid.jl")
include("build_nabla2_matrix.jl")

function test_multiplication()
    Nx = 20
    Ny = 20
    lfgrid = LF2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )

    println("Nx = ", Nx)
    println("Ny = ", Ny)

    ∇2 = build_nabla2_matrix( lfgrid )

    psi1 = rand( lfgrid.Nx, lfgrid.Ny )
    psi2 = rand( lfgrid.Npoints )

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



