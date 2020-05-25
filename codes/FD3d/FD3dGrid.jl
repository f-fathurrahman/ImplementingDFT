include("../FD1d/init_FD1d_grid.jl")
include("../FD1d/init_FD1d_p_grid.jl")

struct FD3dGrid
    Npoints::Int64
    Nx::Int64
    Ny::Int64
    Nz::Int64

    hx::Float64
    hy::Float64
    hz::Float64

    dVol::Float64
    
    x::Array{Float64,1}
    y::Array{Float64,1}
    z::Array{Float64,1}
    
    r::Array{Float64,2}
    
    idx_ip2xyz::Array{Int64,2}
    idx_xyz2ip::Array{Int64,3}

    pbc::Tuple{Bool,Bool,Bool}
end


function FD3dGrid(
    x_domain::Tuple{Float64,Float64}, Nx::Int64,
    y_domain::Tuple{Float64,Float64}, Ny::Int64,
    z_domain::Tuple{Float64,Float64}, Nz::Int64;
    pbc=(false,false,false)
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

    if pbc[3]
        z, hz = init_FD1d_p_grid(z_domain, Nz)        
    else
        z, hz = init_FD1d_grid(z_domain, Nz)
    end

    Npoints = Nx*Ny*Nz
    
    r = zeros(3,Npoints)

    dVol = hx*hy*hz

    idx_ip2xyz = zeros(Int64,3,Npoints)
    idx_xyz2ip = zeros(Int64,Nx,Ny,Nz)

    ip = 0
    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        ip = ip + 1
        r[1,ip] = x[i]
        r[2,ip] = y[j]
        r[3,ip] = z[k]
        idx_ip2xyz[1,ip] = i
        idx_ip2xyz[2,ip] = j
        idx_ip2xyz[3,ip] = k
        idx_xyz2ip[i,j,k] = ip
    end
    
    return FD3dGrid(Npoints, Nx, Ny, Nz, hx, hy, hz, dVol,
        x, y, z, r, idx_ip2xyz, idx_xyz2ip, pbc)
    
end

function FD3dGrid( NN::Array{Int64,1}, AA::Array{Float64,1}, BB::Array{Float64,1}; kwargs... )
    return FD3dGrid( (AA[1], BB[1]), NN[1],
                     (AA[2], BB[2]), NN[2],
                     (AA[3], BB[3]), NN[3]; kwargs... ) 
end