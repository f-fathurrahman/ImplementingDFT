using Printf
using LaTeXStrings

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("init_FD1d_grid.jl")
include("build_D2_matrix_3pt.jl")
include("build_D2_matrix_5pt.jl")
include("build_D2_matrix_7pt.jl")
include("build_D2_matrix_9pt.jl")
include("build_D2_matrix_11pt.jl")

function my_gaussian(x, α=1.0)
    return exp(-α*x^2)
end

function d2_my_gaussian(x, α=1.0)
    return (-2*α + 4*α^2 * x^2) * exp(-α*x^2)
end

function main()
    N = 15
    A = -5.0
    B =  5.0
    x, h = init_FD1d_grid( A, B, N )
    fx = my_gaussian.(x)

    # To emulate continuous function
    Ndense = 200
    x_dense = range(A, stop=B, length=Ndense)
    fx_dense = my_gaussian.(x_dense)
    d2_fx_dense = d2_my_gaussian.(x_dense)

    d2_fx = d2_my_gaussian.(x)

    D2_3pt = build_D2_matrix_3pt(N, h)
    d2_fx_3pt = D2_3pt*fx

    D2_5pt = build_D2_matrix_5pt(N, h)
    d2_fx_5pt = D2_5pt*fx

    D2_7pt = build_D2_matrix_7pt(N, h)
    d2_fx_7pt = D2_7pt*fx

    D2_9pt = build_D2_matrix_9pt(N, h)
    d2_fx_9pt = D2_9pt*fx

    D2_11pt = build_D2_matrix_11pt(N, h)
    d2_fx_11pt = D2_11pt*fx

    plot_title = latexstring("Analytic vs Finite Difference: \$N = $N\$")
    plt.clf()
    #plt.plot(x, d2_fx_3pt, label="FD-3pt", marker="o")
    #plt.plot(x, d2_fx_5pt, label="FD-5pt", marker="o")
    #plt.plot(x, d2_fx_7pt, label="FD-7pt", marker="o")
    #plt.plot(x, d2_fx_9pt, label="FD-9pt", marker="o")
    plt.plot(x, d2_fx_11pt, label="FD-11pt", marker="o")
    plt.plot(x, d2_fx, label="analytic")
    plt.legend()
    plt.title(plot_title)
    plt.tight_layout()
    plt.savefig("IMG_compare_d2_gaussian_"*string(N)*".pdf")


    plot_title = latexstring("Finite Difference - Analytic: \$N = $N\$")
    plt.clf()
    #plt.plot(x, d2_fx_3pt - d2_fx, label="FD-3pt", marker="o")
    #plt.plot(x, d2_fx_5pt - d2_fx, label="FD-5pt", marker="o")
    #plt.plot(x, d2_fx_7pt - d2_fx, label="FD-7pt", marker="o")
    #plt.plot(x, d2_fx_9pt - d2_fx, label="FD-9pt", marker="o")
    plt.plot(x, d2_fx_11pt - d2_fx, label="FD-11pt", marker="o")
    plt.legend()
    plt.title(plot_title)
    plt.tight_layout()
    plt.savefig("IMG_compare_d2_gaussian_error_"*string(N)*".pdf")


    @printf("N = %d\n", N)
    @printf("h = %f\n", h)

    println()

    err = sum(abs.(d2_fx_3pt - d2_fx))/N
    @printf("3pt avg abs err  = %18.10e\n", err)

    err = sum(abs.(d2_fx_5pt - d2_fx))/N
    @printf("5pt avg abs err  = %18.10e\n", err)

    err = sum(abs.(d2_fx_7pt - d2_fx))/N
    @printf("7pt avg abs err  = %18.10e\n", err)

    err = sum(abs.(d2_fx_9pt - d2_fx))/N
    @printf("9pt avg abs err  = %18.10e\n", err)

    err = sum(abs.(d2_fx_11pt - d2_fx))/N
    @printf("11pt avg abs err = %18.10e\n", err)

    println()

    err = sqrt( sum( (d2_fx_3pt - d2_fx).^2 )/N )
    @printf("3pt RMS err  = %18.10e\n", err)

    err = sqrt( sum( (d2_fx_5pt - d2_fx).^2 )/N )
    @printf("5pt RMS err  = %18.10e\n", err)

    err = sqrt( sum( (d2_fx_7pt - d2_fx).^2 )/N )
    @printf("7pt RMS err  = %18.10e\n", err)

    err = sqrt( sum( (d2_fx_9pt - d2_fx).^2 )/N )
    @printf("9pt RMS err  = %18.10e\n", err)

    err = sqrt( sum( (d2_fx_11pt - d2_fx).^2 )/N )
    @printf("11pt RMS err = %18.10e\n", err)

end

main()
