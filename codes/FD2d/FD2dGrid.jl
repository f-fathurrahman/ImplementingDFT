include("../FD1d/init_FD1d_grid.jl")
include("../FD1d/init_FD1d_p_grid.jl")

struct FD2dGrid
    Npoints::Int64
    #
    Lx::Float64
    Ly::Float64
    #
    Nx::Int64
    Ny::Int64
    #
    hx::Float64
    hy::Float64
    dA::Float64
    #
    x::Array{Float64,1}
    y::Array{Float64,1}
    r::Array{Float64,2}
    #
    idx_ip2xy::Array{Int64,2}
    idx_xy2ip::Array{Int64,2}
    #
    pbc::Tuple{Bool,Bool}
end


function FD2dGrid(
    x_domain::Tuple{Float64,Float64}, Nx::Int64,
    y_domain::Tuple{Float64,Float64}, Ny::Int64;
    pbc=(false,false)
)

    if pbc[1]
        x, hx = init_FD1d_p_grid(x_domain, Nx)
    else
        x, hx = init_FD1d_grid(x_domain, Nx)
    end

    if pbc[2]
        y, hy = init_FD1d_p_grid(y_domain, Ny)
    else
        y, hy = init_FD1d_grid(y_domain, Ny)
    end

    Lx = x_domain[2] - y_domain[1]
    Ly = y_domain[2] - y_domain[1]

    dA = hx*hy
    
    Npoints = Nx*Ny
    r = zeros(2,Npoints)
    ip = 0
    idx_ip2xy = zeros(Int64,2,Npoints)
    idx_xy2ip = zeros(Int64,Nx,Ny)
    for j in 1:Ny
        for i in 1:Nx
            ip = ip + 1
            r[1,ip] = x[i]
            r[2,ip] = y[j]
            idx_ip2xy[1,ip] = i
            idx_ip2xy[2,ip] = j
            idx_xy2ip[i,j] = ip
        end
    end
    
    return FD2dGrid(Npoints, Lx, Ly, Nx, Ny, hx, hy, dA, x, y, r, idx_ip2xy, idx_xy2ip, pbc)
    
end

function FD2dGrid( NN::Array{Int64,1}, AA::Array{Float64,1}, BB::Array{Float64,1}; kwargs... )
    return FD2dGrid( (AA[1], BB[1]), NN[1],
                     (AA[2], BB[2]), NN[2]; kwargs... ) 
end

import Base: show
function show( io::IO, grid::FD2dGrid )

    @printf("-----------------\n")
    @printf("FD2dGrid instance\n")
    @printf("-----------------\n")

    @printf(io, "Nx = %8d, hx = %18.10f\n", grid.Nx, grid.hx)
    @printf(io, "Ny = %8d, hy = %18.10f\n", grid.Ny, grid.hy)    
    @printf(io, "dA = %18.10f\n", grid.dA)
end
show( grid::FD2dGrid ) = show(stdout, grid)
