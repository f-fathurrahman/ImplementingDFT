using Printf
using LinearAlgebra
using SparseArrays
using SpecialFunctions

using AlgebraicMultigrid
using Random
using Serialization

include("INC_sch_3d_LF.jl")
include("H_potential.jl")

function main()

    Random.seed!(1234)

    Nx = 50
    Ny = 50
    Nz = 50
    A = 6.0
    grid = LF3dGrid( (-A,A), Nx, (-A,A), Ny, (-A,A), Nz,
                       type_x=:sinc, type_y=:sinc, type_z=:sinc )

    ∇2 = build_nabla2_matrix( grid )

    Vpot = pot_H_atom( grid )
    #Vpot = pot_Hps_HGH( grid )

    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    println("Building preconditioner")
    prec = aspreconditioner(ruge_stuben(Ham))
    println("Done building preconditioner")

    Nstates = 1  # only choose the lowest lying state
    Npoints = Nx*Ny*Nz
    X = rand(Float64, Npoints, Nstates)
    ortho_sqrt!(X)

    evals = diag_LOBPCG!( Ham, X, prec, verbose=true )

    @printf("\n\nEigenvalues\n")
    for i in 1:Nstates
        @printf("%5d %18.10f\n", i, evals[i])
    end
end

main()



