using Printf
using LaTeXStrings

import PyPlot
const plt = PyPlot

plt.rc("text", usetex=true)

include("init_FD1d_p_grid.jl")
include("build_D2_matrix_p_3pt.jl")
include("build_D2_matrix_p_5pt.jl")
include("build_D2_matrix_p_7pt.jl")
include("build_D2_matrix_p_9pt.jl")
include("build_D2_matrix_p_11pt.jl")

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
    x, h = init_FD1d_p_grid( x_min, x_max, N )
    T = (x_max - x_min)
    ω = 2*pi/T
    fx = my_func.(x, ω=ω)

    Ndense = 200
    x_dense = range(x_min, stop=x_max, length=Ndense)
    fx_dense = my_func.(x_dense, ω=ω)
    d2_fx_dense = d2_my_func.(x_dense, ω=ω)
    
    D2 = build_D2_matrix_p_3pt(N, h)
    d2_fx_3pt = D2*fx

    D2 = build_D2_matrix_p_11pt(N, h)
    d2_fx_11pt = D2*fx

    plot_title = latexstring("\$N = $N\$")

    plt.clf()
    plt.plot(x, fx, marker="o", label=L"Sampled $f(x)$")
    plt.plot(x_dense, fx_dense, label=L"f(x)")
    plt.legend()
    plt.grid()
    plt.title(plot_title)
    plt.savefig("IMG_gauss_sin_"*string(N)*".pdf")
    
    plt.clf()
    plt.plot(x, d2_fx_3pt, marker="o", label=L"Approx $f''(x)$ 3pt")
    plt.plot(x, d2_fx_11pt, marker="o", label=L"Approx $f''(x)$ 11pt")
    plt.plot(x_dense, d2_fx_dense, label=L"f''(x)")
    plt.legend()
    plt.grid()
    plt.title(plot_title)
    plt.savefig("IMG_d2_gauss_sin_"*string(N)*".pdf")

    plt.clf()
    plt.plot(x, log.(abs.(d2_fx_3pt - d2_my_func.(x, ω=ω))), marker="o", label=L"Abs error $f''(x)$ 3pt")
    plt.plot(x, log.(abs.(d2_fx_11pt - d2_my_func.(x, ω=ω))), marker="o", label=L"Abs error $f''(x)$ 11pt")
    plt.legend()
    plt.grid()
    plt.title(plot_title)
    plt.savefig("IMG_d2_gauss_sin_"*string(N)*"_error.pdf")

end

main(15)
main(51)
