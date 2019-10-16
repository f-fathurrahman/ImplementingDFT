using Printf
using LinearAlgebra
using PGFPlotsX
using LaTeXStrings

include("init_FD1d_grid.jl")
include("build_D2_matrix_3pt.jl")

function pot_harmonic( x; ω=1.0 )
    return 0.5 * ω^2 * x^2
end

function main()
    # Initialize the grid points
    xmin = -5.0
    xmax =  5.0
    N = 51
    x, h = init_FD1d_grid(xmin, xmax, N)
    # Build 2nd derivative matrix
    D2 = build_D2_matrix_3pt(N, h)
    # Potential
    Vpot = pot_harmonic.(x)
    # Hamiltonian
    Ham = -0.5*D2 + diagm( 0 => Vpot )
    # Solve the eigenproblem
    evals, evecs = eigen( Ham )
    # We will show the 5 lowest eigenvalues
    Nstates = 5
    @printf("Eigenvalues\n")
    for i in 1:Nstates
        @printf("%5d %18.10f\n", i, evals[i])
    end

    # normalize the first three eigenstates
    for i in 1:3
        ss = dot(evecs[:,i], evecs[:,i])*h
        evecs[:,i] = evecs[:,i]/sqrt(ss)
    end

    # Plot up to 3rd eigenstate
    f = @pgf Axis({ title="N="*string(N), height="10cm", width="15cm", xmajorgrids, ymajorgrids },
        PlotInc(Coordinates(x, evecs[:,1])),
        LegendEntry("1st eigenstate"),
        PlotInc(Coordinates(x, evecs[:,2])),
        LegendEntry("2nd eigenstate"),
        PlotInc(Coordinates(x, evecs[:,3])),
        LegendEntry("3rd eigenstate"),
    )
    pgfsave("IMG_main_harmonic_01_"*string(N)*".pdf", f)
end

main()