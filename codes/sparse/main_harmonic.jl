using Printf
using LinearAlgebra
using SparseArrays
using IncompleteLU
using Random

include("../FD2d/FD2dGrid.jl")
include("../FD2d/build_nabla2_matrix.jl")

include("../LF2d/LF2dGrid.jl")
include("../LF2d/build_nabla2_matrix.jl")

include("../common/supporting_functions.jl")
include("../common/ortho_sqrt.jl")
include("../common/ortho_gram_schmidt.jl")

include("ilu0.jl")
include("diag_Emin_PCG.jl")
include("diag_LOBPCG.jl")

function pot_harmonic( grid; ω=1.0 )
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

    Nx = 30
    Ny = 30
    #grid = FD2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )
    grid = LF2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny, types=(:sinc,:sinc) )

    ∇2 = build_nabla2_matrix( grid )

    Vpot = pot_harmonic( grid )
    
    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    # may choose between these two
    #prec = ilu(-0.5*∇2)
    prec = ilu(Ham) # this should result in faster convergence

    #@printf("sizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)
    #@printf("sizeof prec = %18.10f MiB\n", Base.summarysize(prec)/1024/1024)

    dVol = grid.dVol
    Nstates = 10
    Npoints = Nx*Ny
    X = rand(Float64, Npoints, Nstates)
    ortho_sqrt!(X, dVol)

    println("Check normalization")
    for i in 1:Nstates
        @printf("dot: %18.10f\n", dot(X[:,i],X[:,i])*dVol)
    end

    println("Check orthogonal")
    ist = 1
    for i in 1:Nstates
        if i != ist
            @printf("dot: %18.10f\n", dot(X[:,i],X[:,ist])*dVol)
        end
    end

    evals = diag_Emin_PCG!( Ham, X, prec, verbose=true )
    
    #evals = diag_LOBPCG!( Ham, X, prec, verbose=true )
    
    X = X/sqrt(grid.dVol) # renormalize

    println("Check normalization")
    for i in 1:Nstates
        @printf("dot: %18.10f\n", dot(X[:,i],X[:,i])*dVol)
    end

    println("Check orthogonal")
    ist = 1
    for i in 1:Nstates
        if i != ist
            @printf("dot: %18.10f\n", dot(X[:,i],X[:,ist])*dVol)
        end
    end

    @printf("\n\nEigenvalues\n")
    for i in 1:Nstates
        @printf("%5d %18.10f\n", i, evals[i])
    end

end

main()



