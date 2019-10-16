using Printf
using PGFPlotsX
using LaTeXStrings

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
    f = @pgf Axis({ title=plot_title, height="10cm", width="15cm", xmajorgrids, ymajorgrids },
        PlotInc(Coordinates(x, fx)),
        LegendEntry(L"Sampled $f(x)$"),
        PlotInc({mark="none"}, Coordinates(x_dense, fx_dense)),
        LegendEntry(L"f(x)"),
        PlotInc({mark="none"}, Coordinates(x_dense, d2_fx_dense)),
        LegendEntry(L"f''(x)"),
        PlotInc(Coordinates(x, d2_fx_3pt)),
        LegendEntry(L"Approx $f''(x)$"),        
    )
    pgfsave("IMG_gaussian_"*string(N)*".pdf", f)
end

main(15)
main(51)
