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

function main()
    # Initialize the grid points
    xmin = -5.0
    xmax =  5.0
    N = 101
    xgrid, dx = init_FD1d_grid(xmin, xmax, N)
    # Build 2nd derivative matrix
    D2 = build_D2_matrix_11pt(N, dx)
    # Potential
    Vpot = pot_gaussian.(xgrid, A=5.0, x0=-1.0) + pot_gaussian.(xgrid, A=5.0, x0=1.0)
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

    plot_title = "N="*string(N)
    plt.clf()
    plt.plot(xgrid, evecs[:,1], label="1st eigenstate")
    plt.plot(xgrid, evecs[:,2], label="2nd eigenstate")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.xlim(-3.0, 3.0)
    plt.savefig("IMG_main_2gaussian_01_"*string(N)*".png")

    plt.clf()
    plt.plot(xgrid, evecs[:,1].^2, label="Bonding")
    plt.plot(xgrid, evecs[:,2].^2, label="Antibonding")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.xlim(-3.0, 3.0)
    plt.savefig("IMG_main_2gaussian_01_rho_"*string(N)*".png")

end

main()