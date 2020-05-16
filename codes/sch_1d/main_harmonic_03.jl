using Printf
using LinearAlgebra
using LaTeXStrings

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("INC_sch_1d.jl")

function pot_harmonic( x; ω=1.0 )
    return 0.5 * ω^2 * x^2
end

function solve_eigvals(N::Int64, D2_func)
    # Initialize the grid points
    xmin = -5.0
    xmax =  5.0
    x, h = init_FD1d_grid(xmin, xmax, N)
    # Build 2nd derivative matrix
    D2 = D2_func(N, h)
    # Potential
    Vpot = pot_harmonic.(x)
    # Hamiltonian
    Ham = -0.5*D2 + diagm( 0 => Vpot )
    # Solve for the eigenvalues only
    evals = eigvals( Ham )
    
    return evals
end

function main()
    Npoints = 40
    funcs = [build_D2_matrix_3pt, build_D2_matrix_5pt,
             build_D2_matrix_7pt, build_D2_matrix_9pt,
             build_D2_matrix_11pt]

    ω = 1.0
    hbar = 1.0
    ist = 1
    E_ana = (2*ist - 1)*ω*hbar/2

    for f in funcs
        evals = solve_eigvals(Npoints, f)
        dE = abs(evals[ist]-E_ana)
        @printf("%20s: %18.10f %18.10e\n", string(f), evals[ist], dE)
    end
end

main()