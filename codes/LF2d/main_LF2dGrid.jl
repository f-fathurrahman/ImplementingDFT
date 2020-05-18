using Printf

include("LF2dGrid.jl")

function main()
    Nx = 3
    Ny = 4
    lfgrid = LF2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )
    println(lfgrid)
    println(lfgrid.x)
    println(lfgrid.y)
    for ip = 1:lfgrid.Npoints
        @printf("%3d %8.3f %8.3f\n", ip, lfgrid.r[1,ip], lfgrid.r[2,ip])
    end

    println()
    for i = 1:Nx
        for j = 1:Ny
            ip = lfgrid.idx_xy2ip[i,j]
            @printf("[%8.3f,%8.3f] ", lfgrid.r[1,ip], lfgrid.r[2,ip])
        end
        @printf("\n")
    end
end

main()