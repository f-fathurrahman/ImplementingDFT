using Printf

import PyPlot
const plt = PyPlot

include("init_FD1d_grid.jl")
include("build_D2_matrix_3pt.jl")
include("build_D2_matrix_5pt.jl")
include("build_D2_matrix_7pt.jl")
include("build_D2_matrix_9pt.jl")

function my_gaussian(x, α=1.0)
    return exp(-α*x^2)
end

function d2_my_gaussian(x, α=1.0)
    return -2*α*exp(-α*x^2) + 4*α^2 * x^2 * exp(-α*x^2)
end

function main()
    N = 50
    x, h = init_FD1d_grid( (-5.0, 5.0), N )
    fx = my_gaussian.(x)

    plt.clf()
    plt.plot(x, fx, marker="o")
    plt.grid()
    plt.savefig("TEMP_my_gaussian.pdf")

    d2_fx = d2_my_gaussian.(x)

    D2_3pt = build_D2_matrix_3pt(N, h)
    D2_fx_3pt = D2_3pt*fx

    D2_5pt = build_D2_matrix_5pt(N, h)
    D2_fx_5pt = D2_5pt*fx

    D2_7pt = build_D2_matrix_7pt(N, h)
    D2_fx_7pt = D2_7pt*fx

    D2_9pt = build_D2_matrix_9pt(N, h)
    D2_fx_9pt = D2_9pt*fx

    plt.clf()
    plt.plot(x, d2_fx, label="analytic")
    plt.plot(x, D2_fx_3pt, label="FD-3pt")    
    plt.plot(x, D2_fx_5pt, label="FD-5pt")
    plt.plot(x, D2_fx_7pt, label="FD-7pt")
    plt.plot(x, D2_fx_7pt, label="FD-9pt")    
    plt.grid()
    plt.legend()
    plt.savefig("TEMP_d2_my_gaussian.pdf")


    plt.clf()
    plt.plot(x, D2_fx_3pt - d2_fx, label="FD-3pt")  
    plt.plot(x, D2_fx_5pt - d2_fx, label="FD-5pt")
    plt.plot(x, D2_fx_7pt - d2_fx, label="FD-7pt")
    plt.plot(x, D2_fx_9pt - d2_fx, label="FD-9pt")
    plt.grid()
    plt.legend()
    plt.savefig("TEMP_d2_my_gaussian_diff.pdf")

    @printf("N = %d\n", N)
    @printf("h = %f\n", h)

    println()

    err = sum(abs.(D2_fx_3pt - d2_fx))/N
    @printf("3pt avg abs err = %18.10e\n", err)

    err = sum(abs.(D2_fx_5pt - d2_fx))/N
    @printf("5pt avg abs err = %18.10e\n", err)

    err = sum(abs.(D2_fx_7pt - d2_fx))/N
    @printf("7pt avg abs err = %18.10e\n", err)

    err = sum(abs.(D2_fx_9pt - d2_fx))/N
    @printf("9pt avg abs err = %18.10e\n", err)

    println()

    err = sqrt( sum( (D2_fx_3pt - d2_fx).^2 )/N )
    @printf("3pt RMS err = %18.10e\n", err)

    err = sqrt( sum( (D2_fx_5pt - d2_fx).^2 )/N )
    @printf("5pt RMS err = %18.10e\n", err)

    err = sqrt( sum( (D2_fx_7pt - d2_fx).^2 )/N )
    @printf("7pt RMS err = %18.10e\n", err)

    err = sqrt( sum( (D2_fx_9pt - d2_fx).^2 )/N )
    @printf("9pt RMS err = %18.10e\n", err)

end

main()
