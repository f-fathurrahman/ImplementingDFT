using PGFPlotsX
using NLsolve

include("../1d/init_FD1d_grid.jl")
include("deformation_func_sech.jl")

my_deformation_func(r) = deformation_func(r, 2.5, 2.0)

# 1d case
function mapping_func!(F, x, ξ, R)
    Npoints = length(F)
    Ncenter = length(R)
    for i in 1:Npoints
        xf = 0.0
        for α in 1:Ncenter
            r = abs(x[i] - R[α]) # should be changed to distance
            xf = xf + (x[i] - R[α])*my_deformation_func(r)
        end
        F[i] = x[i] + xf - ξ[i]
    end
    return
end


function my_gaussian(x; α=1.0, x0=0.0)
    return exp(-α*(x - x0)^2)
end

function d2_my_gaussian(x; α=1.0, x0=0.0)
    return -2*α*exp(-α*(x-x0)^2) + 4*α^2 * (x-x0)^2 * exp(-α*(x-x0)^2)
end


function main()

    A = -10.0
    B =  10.0
    Npoints = 45
    
    # regularly spaced points
    ξ, h_ξ = init_FD1d_grid(A, B, Npoints)

    centers = [0.0]
    my_mapping_func!(F, x) = mapping_func!(F, x, ξ, centers)

    x = collect(range(-5, stop=5, length=Npoints))
    F = zeros(Npoints)

    initial_x = copy(x)
    sol = nlsolve(my_mapping_func!, initial_x, show_trace=true, ftol=1e-10)

    x_adapt = copy(sol.zero)

    println("x_adapt first = ", x_adapt[1])
    println("x_adapt last  = ", x_adapt[end])

    println(ξ[1])
    scal = x_adapt[1]/ξ[1]
    println(scal)

    y = ones(Npoints)
    y_ξ = -ones(Npoints) 
    y2 = zeros(Npoints)

    f = @pgf Axis({height="10cm", width="15cm", ymin=-3.0, ymax=3.0, xmajorgrids, },
        PlotInc(Coordinates(x_adapt, y)),
        PlotInc(Coordinates(ξ, y_ξ)),
    )
    pgfsave("TEMP_x_adapt.pdf", f)

    c_my_gaussian(x) = my_gaussian(x, x0=0.0, α=10.0)
    c_d2_my_gaussian(x) = d2_my_gaussian(x, x0=0.0, α=10.0)

    fg = c_my_gaussian.(x_adapt)
    d2_fg = c_d2_my_gaussian.(x_adapt)

    x_regular = range(x_adapt[1], stop=x_adapt[end], length=Npoints)
    x_dense = range(x_adapt[1], stop=x_adapt[end], length=20*Npoints)

    fg_regular = c_my_gaussian.(x_regular)
    d2_fg_regular = c_d2_my_gaussian.(x_regular)

    fg_dense = c_my_gaussian.(x_dense)
    d2_fg_dense = c_d2_my_gaussian.(x_dense)

    # plot only at the vicinity of x0
    f = @pgf Axis({height="10cm", width="15cm", xmin=-2, xmax=2, xmajorgrids, ymajorgrids, },
        #
        PlotInc(Coordinates(x_adapt, fg)),
        LegendEntry("fg adapt"),
        #
        PlotInc(Coordinates(x_regular, fg_regular)),
        LegendEntry("fg regular"),
        #
        PlotInc({mark="none"}, Coordinates(x_dense, fg_dense)),
        LegendEntry("fg dense"),
    )
    pgfsave("TEMP_fg_adapt_vs_regular_vs_dense.pdf", f)


    # plot only at the vicinity of x0
    f = @pgf Axis({height="10cm", width="15cm", xmin=-2, xmax=2, xmajorgrids, ymajorgrids, },
        #
        PlotInc(Coordinates(x_adapt, d2_fg)),
        LegendEntry("d2 fg adapt"),
        #
        PlotInc(Coordinates(x_regular, d2_fg_regular)),
        LegendEntry("d2 fg regular"),
        #
        PlotInc({mark="none"}, Coordinates(x_dense, d2_fg_dense)),
        LegendEntry("d2 fg dense"),
    )
    pgfsave("TEMP_d2_fg_adapt_vs_regular_vs_dense.pdf", f)

end

main()