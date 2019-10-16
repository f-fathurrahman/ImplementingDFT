using PGFPlotsX

include("deformation_func.jl")

function main()
    r = range(0.0, stop=8.0, length=100)
    
    fr1 = deformation_func.(r, 2.0, 0.5)
    fr2 = deformation_func.(r, 2.0, 1.5)
    fr3 = deformation_func.(r, 2.0, 2.5)
    fr4 = deformation_func.(r, 2.0, 3.5)

    f = @pgf Axis({title="Varying a"},
        PlotInc({mark="none"}, Coordinates(r,fr1)),
        LegendEntry("a=0.5"),
        PlotInc({mark="none"}, Coordinates(r,fr2)),
        LegendEntry("a=1.5"),
        PlotInc({mark="none"}, Coordinates(r,fr3)),
        LegendEntry("a=2.5"),
        PlotInc({mark="none"}, Coordinates(r,fr4)),
        LegendEntry("a=3.5"),        
    )
    
    pgfsave("TEMP_deformation_func_a.pdf", f)


    fr1 = deformation_func.(r, 1.0, 1.5)
    fr2 = deformation_func.(r, 2.0, 1.5)
    fr3 = deformation_func.(r, 3.0, 1.5)
    fr4 = deformation_func.(r, 4.0, 1.5)

    f = @pgf Axis({title="Varying A"},
        PlotInc({mark="none"}, Coordinates(r,fr1)),
        LegendEntry("A=1.0"),
        PlotInc({mark="none"}, Coordinates(r,fr2)),
        LegendEntry("A=2.0"),
        PlotInc({mark="none"}, Coordinates(r,fr3)),
        LegendEntry("A=3.0"),
        PlotInc({mark="none"}, Coordinates(r,fr4)),
        LegendEntry("A=4.0"),        
    )
    
    pgfsave("TEMP_deformation_func_A.pdf", f)


end

main()