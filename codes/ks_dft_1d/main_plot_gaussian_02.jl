using Printf
using LaTeXStrings

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("INC_sch_1d.jl")

function my_gaussian(x; α=4.0)
    return -10*exp(-α*x^2)
end

function pot(x)
    return my_gaussian(x) + my_gaussian(x - 1.5) + my_gaussian(x + 1.0)
end

function main(N::Int64)
    A = -5.0
    B =  5.0
    xgrid, dx = init_FD1d_grid( A, B, N )
    fx = pot.(xgrid)

    Ndense = 200
    x_dense = range(A, stop=B, length=Ndense)
    fx_dense = pot.(x_dense)

    plot_title = latexstring("\$N = $N\$")
    
    plt.clf()
    plt.plot(xgrid, fx, marker="o", label=L"Sampled $f(x)$")
    plt.plot(x_dense, fx_dense, label=L"f(x)")
    plt.legend()
    plt.grid()
    plt.title(plot_title)
    plt.savefig("IMG_gaussian_f_"*string(N)*".png", dpi=150)
end

main(51)
