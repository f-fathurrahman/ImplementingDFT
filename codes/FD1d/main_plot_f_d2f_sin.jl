using Printf
using LaTeXStrings

import PyPlot
const plt = PyPlot

plt.rc("text", usetex=true)

include("init_FD1d_p_grid.jl")
include("build_D2_matrix_p_3pt.jl")

function my_sin(x; ω=1.0)
    return sin(ω*x)
end

function d2_my_sin(x; ω=1.0)
    return -ω^2*sin(ω*x)
end

function main(N::Int64)
    A = -5.0
    B =  5.0
    x, h = init_FD1d_p_grid( A, B, N )
    ω = 1.5
    fx = my_sin.(x, ω=ω)

    Ndense = 200
    x_dense = range(A, stop=B, length=Ndense)
    fx_dense = my_sin.(x_dense, ω=ω)
    d2_fx_dense = d2_my_sin.(x_dense, ω=ω)
    
    D2 = build_D2_matrix_p_3pt(N, h)
    d2_fx_3pt = D2*fx

    plot_title = latexstring("\$N = $N\$")
    
    plt.clf()
    plt.plot(x, fx, marker="o", label=L"Sampled $f(x)$")
    plt.plot(x_dense, fx_dense, label=L"f(x)")
    plt.plot(x, d2_fx_3pt, marker="o", label=L"Approx $f''(x)$")
    plt.plot(x_dense, d2_fx_dense, label=L"f''(x)")
    plt.legend()
    plt.grid()
    plt.title(plot_title)
    plt.savefig("IMG_sin_"*string(N)*".pdf")

end

main(15)
#main(51)
