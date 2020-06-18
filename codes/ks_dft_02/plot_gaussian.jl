using PGFPlotsX
using LaTeXStrings

function pot_gaussian( r; A=1.0, α=1.0, r0=0.0, normalized=false )
    Npoints = length(r)
    Vpot = zeros(Npoints)
    if normalized
        NN = 1.0/sqrt(α/pi)
    else
        NN = 1.0
    end
    for i in 1:Npoints
        r2 = (r[i] - r0)^2
        Vpot[i] = -A*exp(-α*r2)*sqrt(α/pi)*NN
    end
    return Vpot
end

function main()

    Δ = 0.1

    r = -3.0:Δ:3.0
    Vr1  = pot_gaussian(r)
    Vr2 = pot_gaussian(r, α=3.0)
    Vr3 = pot_gaussian(r, α=5.0)

    println(sum(Vr1)*Δ)
    println(sum(Vr2)*Δ)
    println(sum(Vr3)*Δ)

    fig = @pgf Axis(
        PlotInc(Coordinates(r,Vr1)), LegendEntry(L"\alpha = 1.0"),
        PlotInc(Coordinates(r,Vr2)), LegendEntry(L"\alpha = 3.0"),
        PlotInc(Coordinates(r,Vr3)), LegendEntry(L"\alpha = 5.0"),
    )
    pgfsave("TEMP_gauss.pdf", fig)
end

main()
