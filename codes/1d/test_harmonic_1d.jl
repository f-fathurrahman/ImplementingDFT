using Printf
using LinearAlgebra

include("init_FD1d_grid.jl")
include("build_D2_matrix_3pt.jl")
include("build_D2_matrix_5pt.jl")
include("build_D2_matrix_7pt.jl")
include("build_D2_matrix_9pt.jl")

function pot_harmonic( x; ω=1.0 )
    return 0.5 * ω^2 * x^2
end

function main()

    xmin = -5.0
    xmax =  5.0
    N = 50
    x, h = init_FD1d_grid(xmin, xmax, N)

    D2 = build_D2_matrix_9pt(N, h)
    Vpot = pot_harmonic.(x)
    
    Ham = -0.5*D2 + diagm( 0 => Vpot )

    evals, evecs = eigen( Ham )

    Nstates = 5
    @printf("\n\nEigenvalues\n")
    for i in 1:Nstates
        @printf("%5d %18.10f\n", i, evals[i])
    end
end

main()