using Printf
using LinearAlgebra
using SparseArrays
using IncompleteLU
using Random

import PyPlot
const plt = PyPlot

include("INC_sch_2d_LF.jl")

function pot_harmonic( lfgrid::LF2dGrid; ω=1.0 )
    Npoints = lfgrid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = lfgrid.r[1,i]
        y = lfgrid.r[2,i]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 )
    end
    return Vpot
end

function main()

    Random.seed!(1234)

    Nx = 50
    Ny = 50
    lfgrid = LF2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny, type_x=:sinc, type_y=:sinc )

    ∇2 = build_nabla2_matrix( lfgrid )

    Vpot = pot_harmonic( lfgrid )
    
    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    # may choose between these two
    #prec = ilu(-0.5*∇2)
    prec = ilu(Ham) # this should result in faster convergence

    Nstates = 10
    Npoints = Nx*Ny
    X = rand(Float64, Npoints, Nstates)
    ortho_sqrt!(X)
    evals = diag_Emin_PCG!( Ham, X, prec, verbose=true )
    #evals = diag_LOBPCG!( Ham, X, prec, verbose=true )
    X = X/sqrt(lfgrid.dA) # renormalize

    @printf("\n\nEigenvalues\n")
    for i in 1:Nstates
        @printf("%5d %18.10f\n", i, evals[i])
    end

end

main()



