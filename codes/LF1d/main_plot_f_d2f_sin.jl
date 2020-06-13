using Printf
using LaTeXStrings

import PyPlot
const plt = PyPlot

plt.rc("text", usetex=true)

include("init_LF1d_p_grid.jl")
include("build_D2_matrix_LF1d_p.jl")

function my_sin(x; ω=1.0)
    return sin(ω*x)
end

function d2_my_sin(x; ω=1.0)
    return -ω^2*sin(ω*x)
end

function main(N::Int64)
    x_min = 0.0
    x_max = 2.0
    x, h = init_LF1d_p_grid( x_min, x_max, N )
    L = (x_max - x_min)
    ω = 2*pi/L
    fx = my_sin.(x, ω=ω)

    Ndense = 200
    x_dense = range(x_min, stop=x_max, length=Ndense)
    fx_dense = my_sin.(x_dense, ω=ω)
    d2_fx_dense = d2_my_sin.(x_dense, ω=ω)
    
    D2 = build_D2_matrix_LF1d_p(L, N)    
    d2_fx = D2*fx

    plot_title = latexstring("\$N = $N\$")
    
    plt.clf()
    plt.plot(x, fx, marker="o", label=L"Sampled $f(x)$")
    plt.plot(x_dense, fx_dense, label=L"f(x)")
    plt.plot(x, d2_fx, marker="o", label=L"Approx $f''(x)$")
    plt.plot(x_dense, d2_fx_dense, label=L"f''(x)")
    plt.legend()
    plt.grid()
    plt.title(plot_title)
    plt.savefig("IMG_sin_d2f_"*string(N)*".pdf")

end

main(15)
#main(51)
