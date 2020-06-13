using Printf
using LaTeXStrings

import PyPlot
const plt = PyPlot

plt.rc("text", usetex=true)

include("init_LF1d_p_grid.jl")
include("build_D2_matrix_LF1d_p.jl")

function my_func(x; ω=1.0)
    t = sin(ω*x)
    return exp(-10*t^2) # user larger prefactor to make the peak sharper
end

function d2_my_func(x; ω=1.0)
    return 20*ω^2*(20*sin(ω*x)^2*cos(ω*x)^2 +
        sin(ω*x)^2 - cos(ω*x)^2)*exp(-10*sin(ω*x)^2)
end

function main(N::Int64)
    x_min = 0.0
    x_max = 2.0
    x, h = init_LF1d_p_grid( x_min, x_max, N )
    L = (x_max - x_min)
    ω = 2*pi/L
    fx = my_func.(x, ω=ω)

    Ndense = 200
    x_dense = range(x_min, stop=x_max, length=Ndense)
    fx_dense = my_func.(x_dense, ω=ω)
    d2_fx_dense = d2_my_func.(x_dense, ω=ω)
    
    D2 = build_D2_matrix_LF1d_p(L, N)
    d2_fx = D2*fx

    plot_title = latexstring("\$N = $N\$")

    plt.clf()
    plt.plot(x, fx, marker="o", label=L"Sampled $f(x)$")
    plt.plot(x_dense, fx_dense, label=L"f(x)")
    plt.legend()
    plt.grid()
    plt.title(plot_title)
    plt.savefig("IMG_gauss_sin_"*string(N)*".pdf")
    
    plt.clf()
    plt.plot(x, d2_fx, marker="o", label=L"Approx $f''(x)$")
    plt.plot(x_dense, d2_fx_dense, label=L"f''(x)")
    plt.legend()
    plt.grid()
    plt.title(plot_title)
    plt.savefig("IMG_d2_gauss_sin_"*string(N)*".pdf")

    plt.clf()
    plt.plot(x, log.(abs.(d2_fx - d2_my_func.(x, ω=ω))), marker="o", label=L"log abs error $f''(x)$")
    plt.legend()
    plt.grid()
    plt.title(plot_title)
    plt.savefig("IMG_d2_gauss_sin_"*string(N)*"_error.pdf")

end

main(15)
main(51)
