using Printf

include("FD2dGrid.jl")

function main()
    Nx = 3
    Ny = 4
    fdgrid = FD2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )
    println(fdgrid)
    println(fdgrid.x)
    println(fdgrid.y)
    for ip = 1:fdgrid.Npoints
        @printf("%3d %8.3f %8.3f\n", ip, fdgrid.r[1,ip], fdgrid.r[2,ip])
    end

    println()
    for i = 1:Nx
        for j = 1:Ny
            ip = fdgrid.idx_xy2ip[i,j]
            @printf("[%8.3f,%8.3f] ", fdgrid.r[1,ip], fdgrid.r[2,ip])
        end
        @printf("\n")
    end
end

main()