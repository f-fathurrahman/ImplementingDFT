struct FD1dGrid
    Npoints::Int64
    Lx::Float64
    Nx::Int64
    hx::Float64
    dVol::Float64
    x::Array{Float64,1}
    pbc::Bool
end


function FD1dGrid( x_domain::Tuple{Float64,Float64}, Nx::Int64, pbc=false )

    @assert x_domain[1] < x_domain[2]
    if pbc
        x, hx = init_FD1d_p_grid(x_domain, Nx)
    else
        x, hx = init_FD1d_grid(x_domain, Nx)
    end

    Lx = x_domain[2] - x_domain[1]

    dVol = hx
    
    Npoints = Nx

    return FD1dGrid(Npoints, Lx, Nx, hx, dVol, x, pbc)
    
end

function FD1dGrid( Nx, xi, xf;  kwargs... )
    @assert xi < xf
    return FD1dGrid( (xi, xf), Nx;  kwargs... )
end

import Base: show
function show( io::IO, grid::FD1dGrid )

    @printf("-----------------\n")
    @printf("FD1dGrid instance\n")
    @printf("-----------------\n")

    Nx = grid.Nx; hx = grid.hx

    println()
    @printf(io, "Box size         = %10.5f\n", grid.Lx)
    @printf(io, "Grid spacing     = %10.5f\n", hx)
    @printf(io, "Number of points = %10d\n", grid.Npoints)
    @printf(io, "dVol = dx        = %10.5f\n", grid.dVol)
    println()
    @printf("Some grid points:\n")
    println()
    @printf("%8d %10.5f\n", 1, grid.x[1])
    @printf("%8d %10.5f\n", 2, grid.x[2])
    @printf("       ..   .......\n")
    @printf("%8d %10.5f\n", Nx-1, grid.x[Nx-1])
    @printf("%8d %10.5f\n", Nx, grid.x[Nx])
    println()

    println()
    println("grid.pbc = ", grid.pbc)
    println()

end
show( grid::FD1dGrid ) = show(stdout, grid)
