using Printf
using LinearAlgebra
using PGFPlotsX
using LaTeXStrings

include("init_FD1d_grid.jl")
include("build_D2_matrix_3pt.jl")
include("build_D2_matrix_5pt.jl")
include("build_D2_matrix_7pt.jl")
include("build_D2_matrix_9pt.jl")

# one-dimensional Poeschl-Teller potential
function pot_poeschl_teller( x; λ::Int64=1 )
    return -0.5*λ*(λ + 1)*sech(x)^2
end

function main()
    # Initialize the grid points
    xmin = -5.0
    xmax =  5.0
    N = 51
    x, h = init_FD1d_grid(xmin, xmax, N)
    # Build 2nd derivative matrix
    D2 = build_D2_matrix_9pt(N, h)
    # Potential
    λ = 5
    Vpot = pot_poeschl_teller.(x, λ=λ)
    # Hamiltonian
    Ham = -0.5*D2 + diagm( 0 => Vpot )
    # Solve the eigenproblem
    evals, evecs = eigen( Ham )
    # Show the bound states energies
    Nstates = λ
    @printf("Eigenvalues\n")
    for i in 1:Nstates
        μ = (λ - i + 1)
        E_analytic = -0.5*μ^2
        @printf("%5d %18.10f %18.10f %18.10e\n", i, evals[i], E_analytic, abs(E_analytic-evals[i]))
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
    pgfsave("IMG_main_poeschl_teller_"*string(N)*".pdf", f)
end

main()