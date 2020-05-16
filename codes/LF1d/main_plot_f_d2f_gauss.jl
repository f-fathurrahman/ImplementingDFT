using Printf
using LaTeXStrings

import PyPlot
const plt = PyPlot

plt.rc("text", usetex=true)

include("init_LF1d_c_grid.jl")
include("build_D2_matrix_LF1d_c.jl")

function my_gaussian(x, α=1.0)
    return exp(-α*x^2)
end

function d2_my_gaussian(x, α=1.0)
    return (-2*α + 4*α^2 * x^2) * exp(-α*x^2)
end

function main(N::Int64)
    x_min = -5.0
    x_max =  5.0
    x, h = init_LF1d_c_grid( x_min, x_max, N )
    fx = my_gaussian.(x)

    # Emulate exact 2nd derivative
    Ndense = 200
    x_dense = range(x_min, stop=x_max, length=Ndense)
    fx_dense = my_gaussian.(x_dense)
    d2_fx_dense = d2_my_gaussian.(x_dense)
    
    D2 = build_D2_matrix_LF1d_c(x_min, x_max, N)
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
    plt.savefig("IMG_gaussian_"*string(N)*".pdf")

end

main(15)
main(51)
