# No spherical constraint
struct GVectors
    G::Array{Float64,2}
    G2::Array{Float64,1}
end


function GVectors( grid )

    Npoints = grid.Npoints

    G = zeros(Float64,3,Npoints)
    G2 = zeros(Float64,Npoints)

    Nx = grid.Nx
    Ny = grid.Ny
    Nz = grid.Nz

    Lx = grid.Lx
    Ly = grid.Ly
    Lz = grid.Lz

    ig = 0
    for k in 0:Nz-1, j in 0:Ny-1, i in 0:Nx-1
        #
        ig = ig + 1
        #
        ii = _flip_fft( i, Nx )
        jj = _flip_fft( j, Ny )
        kk = _flip_fft( k, Nz )
        #
        G[1,ig] = ii * 2.0*pi/Lx
        G[2,ig] = jj * 2.0*pi/Ly
        G[3,ig] = kk * 2.0*pi/Lz
        #
        G2[ig] = G[1,ig]^2 + G[2,ig]^2 + G[3,ig]^2
        #
    end

    return GVectors(G, G2)

end

function _flip_fft( mm, S )
    if mm > S/2
        return mm - S
    else
        return mm
    end
end
