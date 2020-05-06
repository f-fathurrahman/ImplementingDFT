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

function main()
    # Initialize the grid points
    xmin = -5.0
    xmax =  5.0
    N = 51
    x, h = init_FD1d_grid(xmin, xmax, N)
    # Build 2nd derivative matrix
    D2 = build_D2_matrix_9pt(N, h)
    # Potential
    Vpot = pot_harmonic.(x)
    # Hamiltonian
    Ham = -0.5*D2 + diagm( 0 => Vpot )
    # Solve the eigenproblem
    evals, evecs = eigen( Ham )
    # We will show the 5 lowest eigenvalues
    Nstates = 5
    @printf("Eigenvalues\n")
    ω = 1.0
    hbar = 1.0
    for i in 1:Nstates
        E_ana = (2*i - 1)*ω*hbar/2
        @printf("%5d %18.10f %18.10f %18.10e\n", i, evals[i], E_ana, abs(evals[i]-E_ana))
    end

end

main()