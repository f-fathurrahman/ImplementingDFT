using Printf
using LaTeXStrings

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("init_LF1d_c_grid.jl")

function my_gaussian(x::Float64; α=1.0)
    return exp( -α*x^2 )
end

function main()
    x_min = -5.0
    x_max =  5.0
    Npoints = 8
    x, h = init_LF1d_c_grid( x_min, x_max, Npoints )
    @printf("Grid spacing = %f\n", h)
    @printf("\nGrid points:\n")
    for i in 1:Npoints
        @printf("%3d %18.10f\n", i, x[i])
    end

    NptsPlot = 200
    x_dense = range(x_min, stop=x_max, length=NptsPlot)

    plt.clf()
    plt.plot(x_dense, my_gaussian.(x_dense), label=L"f(x)")
    plt.plot(x, my_gaussian.(x), label=L"Sampled $f(x)$", marker="o")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig("IMG_gaussian_LF1d_c.pdf")
end

main()
