using Printf
using PGFPlotsX
using LinearAlgebra
using SparseArrays

import PyPlot
const plt = PyPlot

include("FD2dGrid.jl")
include("build_nabla2_matrix.jl")
include("../common/supporting_functions.jl")

function my_gaussian( grid::FD2dGrid; α=1.0 )
    Npoints = grid.Npoints
    f = zeros(Npoints)
    for i in 1:Npoints
        x = grid.r[1,i]
        y = grid.r[2,i]
        r2 = x^2 + y^2
        f[i] = exp(-α*r2)
    end
    return f
end

function main()

    Nx = 75
    Ny = 75
    grid = FD2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )

    ∇2 = build_nabla2_matrix( grid )

    fg = my_gaussian(grid, α=0.5)

    plt.clf()
    plt.surf(grid.x, grid.y, reshape(fg, grid.Nx, grid.Ny), cmap=:jet)
    plt.gca(projection="3d").view_init(30,7)
    fileplot = "IMG_gaussian2d.pdf"
    plt.savefig(fileplot)
    run(`pdfcrop $fileplot $fileplot`)


    d2fg = ∇2*fg

    plt.clf()
    plt.surf(grid.x, grid.y, reshape(d2fg, grid.Nx, grid.Ny), cmap=:jet)
    plt.gca(projection="3d").view_init(30,7)    
    fileplot = "IMG_d2_gaussian2d.pdf"
    plt.savefig(fileplot)
    run(`pdfcrop $fileplot $fileplot`)

end

main()
