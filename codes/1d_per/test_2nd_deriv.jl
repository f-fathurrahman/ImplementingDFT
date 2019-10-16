using Printf

using PGFPlotsX
using LaTeXStrings

include("init_FD1d_grid_per.jl")
include("build_D2_matrix_3pt_per.jl")
include("build_D2_matrix_5pt_per.jl")
#include("build_D2_matrix_7pt.jl")
#include("build_D2_matrix_9pt.jl")

function my_sin(x; period=1.0, α=1.0)
    ω = 2*pi/period
    return sin(ω*x)
end

function d2_my_sin(x; period=1.0, α=1.0)
    ω = 2*pi/period
    return -ω^2*sin(ω*x)
end

function main()
    N = 15
    A = -5.0
    B =  5.0
    L = B - A
    x, h = init_FD1d_grid_per( A, B, N )
    fx = my_sin.(x, period=L)

    Ndense = 200
    x_dense = range(A, stop=B, length=Ndense)
    fx_dense = my_sin.(x_dense, period=L)
    d2_fx_dense = d2_my_sin.(x_dense, period=L)

    plot_title = latexstring("\$N = $N\$")
    f = @pgf Axis({ title=plot_title, height="10cm", width="15cm", xmajorgrids, ymajorgrids },
        PlotInc(Coordinates(x, fx)),
        LegendEntry(L"Sampled $f(x)$"),
        PlotInc({mark="none"}, Coordinates(x_dense, fx_dense)),
        LegendEntry(L"f(x)"),
        PlotInc({mark="none"}, Coordinates(x_dense, d2_fx_dense)),
        LegendEntry(L"f''(x)"),
    )
    pgfsave("TEMP_my_sin_"*string(N)*".pdf", f)

    d2_fx = d2_my_sin.(x, period=L)

    D2_3pt = build_D2_matrix_3pt_per(N, h)
    d2_fx_3pt = D2_3pt*fx

    D2_5pt = build_D2_matrix_5pt_per(N, h)
    d2_fx_5pt = D2_5pt*fx

    #D2_7pt = build_D2_matrix_7pt(N, h)
    #d2_fx_7pt = D2_7pt*fx

    #D2_9pt = build_D2_matrix_9pt(N, h)
    #d2_fx_9pt = D2_9pt*fx

    plot_title = latexstring("Analytic vs Finite Difference: \$N = $N\$")
    f = @pgf Axis({ title=plot_title, height="10cm", width="15cm", xmajorgrids, ymajorgrids },
        PlotInc(Coordinates(x, d2_fx_3pt)),
        LegendEntry("FD-3pt"),
        PlotInc(Coordinates(x, d2_fx_5pt)),
        LegendEntry("FD-5pt"),
        #PlotInc(Coordinates(x, d2_fx_7pt)),
        #LegendEntry("FD-7pt"),
        #PlotInc(Coordinates(x, d2_fx_9pt)),
        #LegendEntry("FD-9pt"),
        PlotInc(Coordinates(x, d2_fx)),
        LegendEntry("analytic"),
    )
    pgfsave("TEMP_d2_my_sin_"*string(N)*".pdf", f)

#=
    plot_title = latexstring("Finite Difference - Analytic: \$N = $N\$")
    f = @pgf Axis({ title=plot_title, height="10cm", width="15cm", xmajorgrids, ymajorgrids },
        PlotInc(Coordinates(x, d2_fx_3pt - d2_fx)),
        LegendEntry("FD-3pt"),
        PlotInc(Coordinates(x, d2_fx_5pt - d2_fx)),
        LegendEntry("FD-5pt"),
        PlotInc(Coordinates(x, d2_fx_7pt - d2_fx)),
        LegendEntry("FD-7pt"),
        PlotInc(Coordinates(x, d2_fx_9pt - d2_fx)),
        LegendEntry("FD-9pt"),
    )
    pgfsave("TEMP_d2_my_gaussian_diff_"*string(N)*".pdf", f)

    @printf("N = %d\n", N)
    @printf("h = %f\n", h)

    println()

    err = sum(abs.(d2_fx_3pt - d2_fx))/N
    @printf("3pt avg abs err = %18.10e\n", err)

    err = sum(abs.(d2_fx_5pt - d2_fx))/N
    @printf("5pt avg abs err = %18.10e\n", err)

    err = sum(abs.(d2_fx_7pt - d2_fx))/N
    @printf("7pt avg abs err = %18.10e\n", err)

    err = sum(abs.(d2_fx_9pt - d2_fx))/N
    @printf("9pt avg abs err = %18.10e\n", err)

    println()

    err = sqrt( sum( (d2_fx_3pt - d2_fx).^2 )/N )
    @printf("3pt RMS err = %18.10e\n", err)

    err = sqrt( sum( (d2_fx_5pt - d2_fx).^2 )/N )
    @printf("5pt RMS err = %18.10e\n", err)

    err = sqrt( sum( (d2_fx_7pt - d2_fx).^2 )/N )
    @printf("7pt RMS err = %18.10e\n", err)

    err = sqrt( sum( (d2_fx_9pt - d2_fx).^2 )/N )
    @printf("9pt RMS err = %18.10e\n", err)
=#

end

main()
