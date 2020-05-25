using Printf
using LinearAlgebra
using SparseArrays
using IncompleteLU
using Random

import PyPlot
const plt = PyPlot

include("INC_sch_2d_LF.jl")

function pot_harmonic( grid::LF2dGrid; ω=1.0 )
    Npoints = grid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = grid.r[1,i]
        y = grid.r[2,i]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 )
    end
    return Vpot
end

function main()

    Random.seed!(1234)

    Nx = 35
    Ny = 35
    grid = LF2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )

    ∇2 = build_nabla2_matrix( grid )

    Vpot = pot_harmonic( grid )
    
    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    # may choose between these two
    #prec = ilu(-0.5*∇2)
    prec = ilu(Ham) # this should result in faster convergence

    Nstates = 10
    Npoints = Nx*Ny
    X = rand(Float64, Npoints, Nstates)
    ortho_sqrt!(X)
    #evals = diag_Emin_PCG!( Ham, X, prec, verbose=true )
    evals = diag_LOBPCG!( Ham, X, prec, verbose=true )
    X = X/sqrt(grid.dA) # renormalize

    @printf("\n\nEigenvalues\n")
    for i in 1:Nstates
        @printf("%5d %18.10f\n", i, evals[i])
    end

    #for i in 1:Nstates
    #    #plt.clf()
    #    #plt.surf(grid.x, grid.y, reshape(X[:,i], grid.Nx, grid.Ny), cmap=:jet)
    #    #plt.tight_layout()
    #    #plt.savefig("IMG_harmonic_psi_"*string(i)*".pdf")
#        ρ = X[:,i].*X[:,i]
#        plt.clf()
#        plt.contourf(grid.x, grid.y, reshape(ρ, grid.Nx, grid.Ny), cmap=:jet)
#        plt.axis("equal")
#        plt.tight_layout()
#        plt.savefig("IMG_harmonic_rho_"*string(i)*".png", dpi=150)
#
#    end

end

main()



