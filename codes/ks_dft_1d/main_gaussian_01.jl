using Printf
using LinearAlgebra
using LaTeXStrings

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("INC_sch_1d.jl")

# assume A and α are positive numbers
function pot_gaussian( x; A=1.0, α=1.0, x0=0.0 )
    return -A*exp( -α*(x-x0)^2 )
end

function main(; do_plot=false)
    # Initialize the grid points
    xmin = -5.0
    xmax =  5.0
    N = 51
    xgrid, dx = init_FD1d_grid(xmin, xmax, N)
    # Build 2nd derivative matrix
    D2 = build_D2_matrix_11pt(N, dx)
    # Potential
    Vpot = pot_gaussian.(xgrid, A=10.0)
    # Hamiltonian
    Ham = -0.5*D2 + diagm( 0 => Vpot )
    # Solve the eigenproblem
    evals, evecs = eigen( Ham )
    # We will show the 5 lowest eigenvalues
    Nstates = 5
    @printf(" State         Energy Levels\n")
    for i in 1:Nstates
        @printf("%5d %18.10f\n", i, evals[i])
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
    if do_plot
        plot_title = "N="*string(N)
        plt.plot(xgrid, evecs[:,1], label="1st eigenstate", marker="o")
        plt.plot(xgrid, evecs[:,2], label="2nd eigenstate", marker="o")
        plt.plot(xgrid, evecs[:,3], label="3rd eigenstate", marker="o")
        plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.savefig("IMG_main_gaussian_01_"*string(N)*".pdf")
    end
end

main(do_plot=false)