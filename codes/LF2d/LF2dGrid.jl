include("../LF1d/init_LF1d_c_grid.jl")

struct LF2dGrid
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
end

# FIXME: Add different BC
function LF2dGrid(
    x_domain::Tuple{Float64,Float64},
    Nx::Int64,
    y_domain::Tuple{Float64,Float64},
    Ny::Int64
)

    x, hx = init_LF1d_c_grid(x_domain, Nx)
    y, hy = init_LF1d_c_grid(y_domain, Ny)

    Lx = x_domain[2] - x_domain[1]
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
    
    return LF2dGrid(Npoints, Lx, Ly, Nx, Ny, hx, hy, dA, x, y, r, idx_ip2xy, idx_xy2ip)
    
end

function LF2dGrid( NN::Array{Int64,1}, AA::Array{Float64,1}, BB::Array{Float64,1} )
    return LF2dGrid( (AA[1], BB[1]), NN[1],
                     (AA[2], BB[2]), NN[2] ) 
end

import Base: show
function show( io::IO, lfgrid::LF2dGrid )

    @printf("-----------------\n")
    @printf("LF2dGrid instance\n")
    @printf("-----------------\n")

    @printf(io, "Nx = %8d, hx = %18.10f\n", lfgrid.Nx, lfgrid.hx)
    @printf(io, "Ny = %8d, hy = %18.10f\n", lfgrid.Ny, lfgrid.hy)    
    @printf(io, "dA = %18.10f\n", lfgrid.dA)
end
show( fdgrid::LF2dGrid ) = show(stdout, lfgrid)
