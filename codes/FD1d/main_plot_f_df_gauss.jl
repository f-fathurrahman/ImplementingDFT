using Printf
using LaTeXStrings

import PyPlot
const plt = PyPlot

plt.rc("text", usetex=true)

include("init_FD1d_grid.jl")
include("build_D1_matrix_3pt.jl")
include("build_D1_matrix_5pt.jl")
include("build_D1_matrix_7pt.jl")
include("build_D1_matrix_9pt.jl")
include("build_D1_matrix_11pt.jl")

function my_gaussian(x; α=1.0)
    return exp(-α*x^2)
end

function d1_my_gaussian(x; α=1.0)
    return -2*α*x*exp(-α*x^2)
end

function main(N::Int64)
    A = -5.0
    B =  5.0
    x, h = init_FD1d_grid( A, B, N )
    fx = my_gaussian.(x)

    Ndense = 200
    x_dense = range(A, stop=B, length=Ndense)
    fx_dense = my_gaussian.(x_dense)
    d1_fx_dense = d1_my_gaussian.(x_dense)
    
    #D1 = build_D1_matrix_3pt(N, h)
    #D1 = build_D1_matrix_5pt(N, h)
    #D1 = build_D1_matrix_7pt(N, h)
    #D1 = build_D1_matrix_9pt(N, h)
    D1 = build_D1_matrix_11pt(N, h)
    d1_fx = D1*fx

    plot_title = latexstring("\$N = $N\$")
    
    plt.clf()
    plt.plot(x, fx, marker="o", label=L"Sampled $f(x)$")
    plt.plot(x_dense, fx_dense, label=L"f(x)")
    plt.plot(x, d1_fx, marker="o", label=L"Approx $f'(x)$")
    plt.plot(x_dense, d1_fx_dense, label=L"f'(x)")
    plt.legend()
    plt.grid()
    plt.title(plot_title)
    plt.savefig("IMG_gaussian_df_"*string(N)*".pdf")

end

main(15)
main(51)
