using Printf
using PGFPlotsX
using LinearAlgebra
using SparseArrays

import PyPlot
const plt = PyPlot

include("FD2dGrid.jl")
include("build_nabla2_matrix.jl")
include("../common/supporting_functions.jl")

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

    fg = my_gaussian(fdgrid, α=0.5)

    plt.clf()
    plt.surf(fdgrid.x, fdgrid.y, reshape(fg, fdgrid.Nx, fdgrid.Ny), cmap=:jet)
    plt.gca(projection="3d").view_init(30,7)
    fileplot = "IMG_gaussian2d.pdf"
    plt.savefig(fileplot)
    run(`pdfcrop $fileplot $fileplot`)


    d2fg = ∇2*fg

    plt.clf()
    plt.surf(fdgrid.x, fdgrid.y, reshape(d2fg, fdgrid.Nx, fdgrid.Ny), cmap=:jet)
    plt.gca(projection="3d").view_init(30,7)    
    fileplot = "IMG_d2_gaussian2d.pdf"
    plt.savefig(fileplot)
    run(`pdfcrop $fileplot $fileplot`)

end

main()
