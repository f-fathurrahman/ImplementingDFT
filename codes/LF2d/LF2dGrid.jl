include("../LF1d/init_LF1d_c_grid.jl")
include("../LF1d/init_LF1d_sinc_grid.jl")

struct LF2dGrid
    Npoints::Int64
    types::Tuple{Symbol,Symbol}
    #
    Lx::Float64
    Ly::Float64
    #
    Nx::Int64
    Ny::Int64
    #
    hx::Float64
    hy::Float64
    dVol::Float64
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

# FIXME: Add different BC
function LF2dGrid(
    x_domain::Tuple{Float64,Float64}, Nx::Int64,
    y_domain::Tuple{Float64,Float64}, Ny::Int64;
    types=(:C,:C)
)

    if types[1] == :C
        x, hx = init_LF1d_c_grid(x_domain, Nx)
    elseif types[2] == :sinc
        x, hx = init_LF1d_sinc_grid(x_domain, Nx)
    else
        error("Unsupported types[1] = ", types[1])
    end

    if types[2] == :C
        y, hy = init_LF1d_c_grid(y_domain, Ny)
    elseif types[2] == :sinc
        y, hy = init_LF1d_sinc_grid(y_domain, Ny)
    else
        error("Unsupported types[2] = ", types[2])
    end

    Lx = x_domain[2] - x_domain[1]
    Ly = y_domain[2] - y_domain[1]

    dVol = hx*hy
    
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
    
    pbc1 = (types[1] == :P)
    pbc2 = (types[2] == :P)

    return LF2dGrid( Npoints, types, Lx, Ly, Nx, Ny, hx, hy, dVol,
        x, y, r, idx_ip2xy, idx_xy2ip, (pbc1,pbc2) )
    
end

function LF2dGrid( NN::Array{Int64,1}, AA::Array{Float64,1}, BB::Array{Float64,1}; kwargs... )
    return LF2dGrid( (AA[1], BB[1]), NN[1],
                     (AA[2], BB[2]), NN[2]; kwargs... ) 
end

import Base: show
function show( io::IO, grid::LF2dGrid )

    @printf("-----------------\n")
    @printf("LF2dGrid instance\n")
    @printf("-----------------\n")

    @printf(io, "Nx = %8d, hx = %18.10f\n", grid.Nx, grid.hx)
    @printf(io, "Ny = %8d, hy = %18.10f\n", grid.Ny, grid.hy)
end
show( fdgrid::LF2dGrid ) = show(stdout, grid)
