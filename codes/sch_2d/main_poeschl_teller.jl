using Printf
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using IncompleteLU
using Random

import PyPlot
const plt = PyPlot

include("INC_sch_2d.jl")

function pot_poeschl_teller( grid::FD2dGrid )
    Npoints = grid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = grid.r[1,i]
        y = grid.r[2,i]
        r = sqrt(x^2 + y^2)
        Vpot[i] = -1/cosh(r)^2
    end
    return Vpot
end

function main()

    Random.seed!(1234)

    Nx = 50
    Ny = 50
    grid = FD2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )

    ∇2 = build_nabla2_matrix( grid )

    Vpot = pot_poeschl_teller( grid )
    
    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    # may choose between these two
    #prec = ilu(-0.5*∇2)
    prec = ilu(Ham) # this should result in faster convergence

    Nstates = 5
    Npoints = Nx*Ny
    X = rand(Float64, Npoints, Nstates)
    ortho_sqrt!(X)
    evals = diag_Emin_PCG!( Ham, X, prec, verbose=true )

    @printf("\n\nEigenvalues\n")
    for i in 1:Nstates
        @printf("%5d %18.10f\n", i, evals[i])
    end
end

main()



