using Printf
using LinearAlgebra
using SparseArrays
using IncompleteLU
using PGFPlotsX
using Random

include("FD2dGrid.jl")
include("build_nabla2_matrix.jl")
include("../supporting_functions.jl")
include("../diag_Emin_PCG.jl")
include("../diag_LOBPCG.jl")
include("../ortho_sqrt.jl")
include("../ortho_gram_schmidt.jl")

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

    ∇2 = build_nabla2_matrix( fdgrid, func_1d=build_D2_matrix_9pt )

    Vpot = pot_harmonic( fdgrid )
    
    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    # may choose between these two
    #prec = ilu(-0.5*∇2)
    prec = ilu(Ham) # this should result in faster convergence

    Nstates = 5
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

    for i in 1:Nstates
        fig = @pgf Axis({ height = "10cm", width = "10cm", view=(20,10), "colormap/jet", },
            Plot3( { surf, },
                Coordinates(fdgrid.x, fdgrid.y, reshape(X[:,i], fdgrid.Nx, fdgrid.Ny) )
            )
        )
        pgfsave("IMG_harmonic_psi_"*string(i)*".pdf", fig)
    end

end

main()



