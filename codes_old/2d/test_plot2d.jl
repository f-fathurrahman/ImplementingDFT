using Printf
using PGFPlotsX
using LinearAlgebra
using SparseArrays

include("FD2dGrid.jl")
include("build_nabla2_matrix.jl")
include("../supporting_functions.jl")

function my_gaussian( fdgrid::FD2dGrid; α=1.0 )
    Npoints = fdgrid.Npoints
    f = zeros(Npoints)
    for i in 1:Npoints
        x = fdgrid.r[1,i]
        y = fdgrid.r[2,i]
        r2 = x^2 + y^2
        f[i] = exp(-α*r2)
    end
    return f
end

function main()

    Nx = 75
    Ny = 75
    fdgrid = FD2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )

    ∇2 = build_nabla2_matrix( fdgrid, func_1d=build_D2_matrix_9pt )

    fg = my_gaussian(fdgrid)

    fig = @pgf Axis({ height = "10cm", width = "10cm", view=(20,10), "colormap/jet" },
        Plot3( { surf },
            Coordinates(fdgrid.x, fdgrid.y, reshape(fg, fdgrid.Nx, fdgrid.Ny) )
        )
    )
    pgfsave("IMG_gaussian2d.pdf", fig)

    d2fg = ∇2*fg
    fig = @pgf Axis({ height = "10cm", width = "10cm", view=(20,10), "colormap/jet", },
        Plot3( { surf, },
            Coordinates(fdgrid.x, fdgrid.y, reshape(d2fg, fdgrid.Nx, fdgrid.Ny) )
        )
    )
    pgfsave("IMG_d2_gaussian2d.pdf", fig)

end

main()
