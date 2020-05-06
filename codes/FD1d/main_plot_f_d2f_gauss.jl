using Printf
using LaTeXStrings

import PyPlot
const plt = PyPlot

plt.rc("text", usetex=true)

include("init_FD1d_grid.jl")
include("build_D2_matrix_3pt.jl")

function my_gaussian(x, α=1.0)
    return exp(-α*x^2)
end

function d2_my_gaussian(x, α=1.0)
    return (-2*α + 4*α^2 * x^2) * exp(-α*x^2)
end

function main(N::Int64)
    A = -5.0
    B =  5.0
    x, h = init_FD1d_grid( A, B, N )
    fx = my_gaussian.(x)

    Ndense = 200
    x_dense = range(A, stop=B, length=Ndense)
    fx_dense = my_gaussian.(x_dense)
    d2_fx_dense = d2_my_gaussian.(x_dense)
    
    D2 = build_D2_matrix_3pt(N, h)
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
    plt.savefig("IMG_gaussian_"*string(N)*".pdf")

end

main(15)
main(51)
