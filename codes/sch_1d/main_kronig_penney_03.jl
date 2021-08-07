using Printf
using LinearAlgebra
using LaTeXStrings

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("INC_sch_1d.jl")

function pot_kronig_penney( x; V0 = 1.0, L=1.0 )
    if x >= L/2
        return V0
    else
        return 0.0
    end
end

function build_Ham_matrix( D1, D2, Vpot, k::Float64 )
    N = size(Vpot,1)
    Ham = -0.5*( D2 + 2*im*k*D1 - k^2*diagm(0 => ones(N)) ) + diagm( 0 => Vpot )
    return Hermitian(Ham)
end

function main()
    # Initialize the grid points
    xmin = 0.0
    xmax = 10.0
    L = xmax - xmin
    N = 101
    x, h = init_FD1d_p_grid(xmin, xmax, N)
    # Build derivative matrices
    D1 = build_D1_matrix_p_11pt(N, h)
    D2 = build_D2_matrix_p_11pt(N, h)
    # Potential
    Vpot = pot_kronig_penney.(x, L=L, V0=0.1)

    Nk = 51
    k = range(-pi/L, pi/L, length=Nk)
    Nstates = 3
    ebands = zeros(Nk, Nstates)
    for ik in 1:Nk
        # Hamiltonian
        Ham = build_Ham_matrix(D1, D2, Vpot, k[ik])
        # Solve the eigenproblem
        evals = eigvals( Ham )
        for ist in 1:Nstates
            ebands[ik,ist] = evals[ist]
        end
    end

    plot_title = "N="*string(N)
    plt.clf()
    plt.figure(figsize=(5,8))
    for ist in 1:Nstates
        labelstr = "band-"*string(ist)
        plt.plot(k, ebands[:,ist], label=labelstr, marker="o")
    end
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig("IMG_main_kronig_03_"*string(N)*".pdf")

    println("At k = ", k[1])
    println(ebands[1,:])

    println("At k = ", k[26])
    println(ebands[26,:])
end

main()