using Printf
using LinearAlgebra
using SparseArrays
using IncompleteLU
using Random

import PyPlot
const plt = PyPlot

include("INC_sch_2d.jl")

function pot_harmonic( fdgrid::FD2dGrid; ω=1.0 )
    Npoints = fdgrid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = fdgrid.r[1,i]
        y = fdgrid.r[2,i]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 )
    end
    return Vpot
end

function main()

    Random.seed!(1234)

    Nx = 50
    Ny = 50
    fdgrid = FD2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )

    ∇2 = build_nabla2_matrix( fdgrid, func_1d=build_D2_matrix_11pt )

    Vpot = pot_harmonic( fdgrid )
    
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
    X = X/sqrt(fdgrid.dA) # renormalize

    @printf("\n\nEigenvalues\n")
    for i in 1:Nstates
        @printf("%5d %18.10f\n", i, evals[i])
    end

    #for i in 1:Nstates
    #    #plt.clf()
    #    #plt.surf(fdgrid.x, fdgrid.y, reshape(X[:,i], fdgrid.Nx, fdgrid.Ny), cmap=:jet)
    #    #plt.tight_layout()
    #    #plt.savefig("IMG_harmonic_psi_"*string(i)*".pdf")
#        ρ = X[:,i].*X[:,i]
#        plt.clf()
#        plt.contourf(fdgrid.x, fdgrid.y, reshape(ρ, fdgrid.Nx, fdgrid.Ny), cmap=:jet)
#        plt.axis("equal")
#        plt.tight_layout()
#        plt.savefig("IMG_harmonic_rho_"*string(i)*".png", dpi=150)
#
#    end

end

main()



