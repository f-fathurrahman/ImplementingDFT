using PGFPlotsX
using LaTeXStrings

include("init_FD1d_grid.jl")

function my_gaussian(x::Float64; α=1.0)
    return exp( -α*x^2 )
end

function main()
    A = -5.0
    B =  5.0
    Npoints = 8
    x, h = init_FD1d_grid( A, B, Npoints )
    println(x)
    println(h)

    NptsPlot = 200
    x_dense = range(A, stop=5, length=NptsPlot)

    f = @pgf( Axis( {height = "6cm", width = "10cm" },
        PlotInc( {mark="none"}, Coordinates(x_dense, my_gaussian.(x_dense)) ),
        LegendEntry(L"f(x)"),
        PlotInc( Coordinates(x, my_gaussian.(x)) ),
        LegendEntry(L"Sampled $f(x)$"),
        )
    )
    pgfsave("TEMP_gaussian_1d.pdf", f)
end

main()