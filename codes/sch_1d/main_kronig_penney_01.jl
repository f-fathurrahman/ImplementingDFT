using Printf
using LinearAlgebra
using LaTeXStrings

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("INC_sch_1d.jl")

function pot_kronig_penney( x; V0 = 1.0, L=1.0 )
    if (x >= L/4) && (x < 3*L/4)
        return 0.0
    else
        return V0
    end
end

function pot_mathieu(x; V0=1.0, L=1.0 )
    return V0*( 1 + cos(2*pi*x/L) )
end

function build_Ham_matrix( D1, D2, Vpot, k::Float64 )
    N = size(Vpot,1)
    Ham = -0.5*( D2 + 2*im*k*D1 - k^2*diagm(0 => ones(N)) ) + diagm( 0 => Vpot )
    return Hermitian(Ham)
end

function main()
    # Initialize the grid points
    xmin = 0.0
    xmax = 5.0
    L = xmax - xmin
    N = 51
    x, h = init_FD1d_p_grid(xmin, xmax, N)
    # Build derivative matrices
    D1 = build_D1_matrix_p_3pt(N, h)
    D2 = build_D2_matrix_p_3pt(N, h)
    # Potential
    Vpot = pot_kronig_penney.(x, L=L, V0=1.5)
    #Vpot = pot_mathieu.(x, L=L, V0=1.5)
    # Hamiltonian
    k = 0.0
    Ham = build_Ham_matrix(D1, D2, Vpot, k)

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
    plot_title = "N="*string(N)
    plt.clf()
    plt.plot(x, real(evecs[:,1]), label="1st eigenstate", marker="o")
    plt.plot(x, real(evecs[:,2]), label="2nd eigenstate", marker="o")
    plt.plot(x, real(evecs[:,3]), label="3rd eigenstate", marker="o")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig("IMG_main_kronig_01_"*string(N)*".pdf")
end

main()