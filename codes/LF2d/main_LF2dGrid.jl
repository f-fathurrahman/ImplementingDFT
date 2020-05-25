using Printf

include("LF2dGrid.jl")

function main_c()
    Nx = 3
    Ny = 4
    grid = LF2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )
    println(grid)
    println(grid.x)
    println(grid.y)
    for ip = 1:grid.Npoints
        @printf("%3d %8.3f %8.3f\n", ip, grid.r[1,ip], grid.r[2,ip])
    end

    println()
    for i = 1:Nx
        for j = 1:Ny
            ip = grid.idx_xy2ip[i,j]
            @printf("[%8.3f,%8.3f] ", grid.r[1,ip], grid.r[2,ip])
        end
        @printf("\n")
    end
end
#main_c()

function main_sinc()
    Nx = 3
    Ny = 4
    grid = LF2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny, types=(:sinc,:sinc) )
    println(grid)
    println(grid.x)
    println(grid.y)
    for ip = 1:grid.Npoints
        @printf("%3d %8.3f %8.3f\n", ip, grid.r[1,ip], grid.r[2,ip])
    end

    println()
    for i = 1:Nx
        for j = 1:Ny
            ip = grid.idx_xy2ip[i,j]
            @printf("[%8.3f,%8.3f] ", grid.r[1,ip], grid.r[2,ip])
        end
        @printf("\n")
    end
end
main_sinc()