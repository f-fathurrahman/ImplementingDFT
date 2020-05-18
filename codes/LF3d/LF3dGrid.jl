include("../LF1d/init_LF1d_c_grid.jl")
include("../LF1d/init_LF1d_sinc_grid.jl")

struct LF3dGrid
    Npoints::Int64
    type_x::Symbol
    type_y::Symbol
    type_z::Symbol

    Lx::Float64
    Ly::Float64
    Lz::Float64

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
end


function LF3dGrid(
    x_domain::Tuple{Float64,Float64}, Nx::Int64,
    y_domain::Tuple{Float64,Float64}, Ny::Int64,
    z_domain::Tuple{Float64,Float64}, Nz::Int64;
    type_x=:C, type_y=:C, type_z=:C
)

    if type_x == :C
        x, hx = init_LF1d_c_grid(x_domain, Nx)
    elseif type_x == :sinc
        x, hx = init_LF1d_sinc_grid(x_domain, Nx)
    else
        error("Unsupported type_x = ", type_x)
    end

    if type_y == :C
        y, hy = init_LF1d_c_grid(y_domain, Ny)
    elseif type_x == :sinc
        y, hy = init_LF1d_sinc_grid(y_domain, Ny)
    else
        error("Unsupported type_y = ", type_y)
    end

    if type_z == :C
        z, hz = init_LF1d_c_grid(z_domain, Nz)
    elseif type_x == :sinc
        z, hz = init_LF1d_sinc_grid(z_domain, Nz)
    else
        error("Unsupported type_z = ", type_z)
    end

    Lx = x_domain[2] - x_domain[1]
    Ly = y_domain[2] - y_domain[1]
    Lz = z_domain[2] - z_domain[1]

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
    
    return LF3dGrid( Npoints, type_x, type_y, type_z, 
        Lx, Lx, Lz, Nx, Ny, Nz, hx, hy, hz, dVol,
        x, y, z, r, idx_ip2xyz, idx_xyz2ip )
    
end

function LF3dGrid( NN::Array{Int64,1}, AA::Array{Float64,1}, BB::Array{Float64,1} )
    return LF3dGrid( (AA[1], BB[1]), NN[1],
                     (AA[2], BB[2]), NN[2],
                     (AA[3], BB[3]), NN[3] ) 
end