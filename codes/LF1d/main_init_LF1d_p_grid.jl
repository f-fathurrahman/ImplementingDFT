using Printf
using LaTeXStrings

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("init_LF1d_p_grid.jl")

function my_sin(x; ω=1.0)
    return sin(ω*x)
end

function main()
    x_min = 0.0
    x_max = 5.0
    Npoints = 9
    x, h = init_LF1d_p_grid( x_min, x_max, Npoints )
    
    T = (x_max - x_min)
    ω = 2*pi/T
    fx = my_sin.(x, ω=ω)
    
    @printf("Grid spacing = %f\n", h)
    @printf("\nGrid points:\n")
    for i in 1:Npoints
        @printf("%3d %18.10f\n", i, x[i])
    end

    Ndense = 200
    x_dense = range(x_min, stop=x_max, length=Ndense)
    fx_dense = my_sin.(x_dense, ω=ω)

    plt.clf()
    plt.plot(x_dense, fx_dense, label=L"f(x)")
    plt.plot(x, fx, label=L"Sampled $f(x)$", marker="o")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig("IMG_sin_LF1d_p.pdf")
end

main()
