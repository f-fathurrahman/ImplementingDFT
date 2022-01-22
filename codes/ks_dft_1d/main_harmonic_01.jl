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
    x, dx = init_FD1d_grid(xmin, xmax, N)
    # Build 2nd derivative matrix
    D2 = build_D2_matrix_11pt(N, dx)
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
    @printf(" State         Approx              Exact          Difference\n")
    for i in 1:Nstates
        E_ana = (2*i - 1)*ω*hbar/2
        @printf("%5d %18.10f %18.10f %18.10e\n", i, evals[i], E_ana, abs(evals[i]-E_ana))
    end

    # normalize the first three eigenstates
    for i in 1:3
        evecs[:,i] = evecs[:,i]/sqrt(dx)
    end

    println("Check ortho:")
    for i in 1:3
        @printf("[%3d,%3d]: %18.10f\n", 1, i, dot(evecs[:,1], evecs[:,i])*dx)
    end

    # Plot up to 3rd eigenstate
    plot_title = "N="*string(N)
    plt.plot(x, evecs[:,1], label="1st eigenstate", marker="o")
    plt.plot(x, evecs[:,2], label="2nd eigenstate", marker="o")
    plt.plot(x, evecs[:,3], label="3rd eigenstate", marker="o")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig("IMG_main_harmonic_01_"*string(N)*".pdf")
end

main()