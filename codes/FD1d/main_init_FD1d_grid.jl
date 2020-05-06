using Printf
using LaTeXStrings

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("init_FD1d_grid.jl")

function my_gaussian(x::Float64; α=1.0)
    return exp( -α*x^2 )
end

function main()
    A = -5.0
    B =  5.0
    Npoints = 8
    x, h = init_FD1d_grid( A, B, Npoints )
    @printf("Grid spacing = %f\n", h)
    @printf("\nGrid points:\n")
    for i in 1:Npoints
        @printf("%3d %18.10f\n", i, x[i])
    end

    NptsPlot = 200
    x_dense = range(A, stop=5, length=NptsPlot)

    plt.clf()
    plt.plot(x_dense, my_gaussian.(x_dense), label=L"f(x)")
    plt.plot(x, my_gaussian.(x), label=L"Sampled $f(x)$", marker="o")
    plt.legend()
    plt.tight_layout()
    plt.savefig("IMG_gaussian_1d.pdf")
end

main()