"""
The type for set of G-vectors for describing density and potentials
"""
struct GVectors
    Ng::Int64
    G::Array{Float64,2}
    G2::Array{Float64,1}
    idx_g2r::Array{Int64,1}
    G2_shells::Array{Float64,1}
    idx_g2shells::Array{Int64,1}
end


function GVectors( grid )

    Npoints = grid.Npoints
    #
    Nx = grid.Nx
    Ny = grid.Ny
    Nz = grid.Nz
    #
    Lx = grid.Lx
    Ly = grid.Ly
    Lz = grid.Lz

    RecVecs = zeros(3,3)
    RecVecs[1,1] = 2.0*pi/Lx
    RecVecs[2,2] = 2.0*pi/Ly
    RecVecs[3,3] = 2.0*pi/Lz

    Δ = max(grid.hx, grid.hy, grid.hz)
    ecutrho = (pi/Δ)^2
    Ns = (Nx,Ny,Nz)
    Ng = calc_Ng( Ns, RecVecs, ecutrho )

    G  = zeros(Float64,3,Ng)
    G2 = zeros(Float64,Ng)
    idx_g2r = zeros(Int64,Ng)

    ig = 0
    ip = 0
    for k in 0:Ns[3]-1, j in 0:Ns[2]-1, i in 0:Ns[1]-1
        ip = ip + 1
        gi = _flip_fft( i, Ns[1] )
        gj = _flip_fft( j, Ns[2] )
        gk = _flip_fft( k, Ns[3] )
        Gx = RecVecs[1,1]*gi + RecVecs[1,2]*gj + RecVecs[1,3]*gk
        Gy = RecVecs[2,1]*gi + RecVecs[2,2]*gj + RecVecs[2,3]*gk
        Gz = RecVecs[3,1]*gi + RecVecs[3,2]*gj + RecVecs[3,3]*gk
        G2_temp = Gx^2 + Gy^2 + Gz^2
        if 0.5*G2_temp <= ecutrho
            ig = ig + 1
            G[1,ig] = Gx
            G[2,ig] = Gy
            G[3,ig] = Gz
            G2[ig] = G2_temp
            idx_g2r[ig] = ip
        end
    end

    # if sorted
    idx_sorted = sortperm(G2)
    G = G[:,idx_sorted]
    G2 = G2[idx_sorted]
    idx_g2r = idx_g2r[idx_sorted]

    G2_shells, idx_g2shells = init_Gshells( G2 )

    return GVectors( Ng, G, G2, idx_g2r, G2_shells, idx_g2shells )

end

"""
Calculates number of G-vectors satisfying |G|^2 <= 2*ecutrho.
This function is used by function `init_gvec`.
"""
function calc_Ng( Ns, RecVecs, ecutrho )
    ig = 0
    Ng = 0
    #
    G = zeros(Float64,3)
    #
    for k in 0:Ns[3]-1, j in 0:Ns[2]-1, i in 0:Ns[1]-1
        ig = ig + 1
        gi = _flip_fft( i, Ns[1] )
        gj = _flip_fft( j, Ns[2] )
        gk = _flip_fft( k, Ns[3] )
        G[1] = RecVecs[1,1]*gi + RecVecs[1,2]*gj + RecVecs[1,3]*gk
        G[2] = RecVecs[2,1]*gi + RecVecs[2,2]*gj + RecVecs[2,3]*gk
        G[3] = RecVecs[3,1]*gi + RecVecs[3,2]*gj + RecVecs[3,3]*gk
        G2 = G[1]^2 + G[2]^2 + G[3]^2
        if 0.5*G2 <= ecutrho
            Ng = Ng + 1
        end
    end
    return Ng
end

function init_Gshells( G2_sorted::Array{Float64,1} )

    eps8 = 1e-8

    Ng = length(G2_sorted)

    ngl = 1
    for ig = 2:Ng
        if G2_sorted[ig] > (G2_sorted[ig-1] + eps8)
            ngl = ngl + 1
        end
    end

    G2_shells = zeros(ngl)
    idx_g2shells = zeros(Int64,Ng)
    
    G2_shells[1] = G2_sorted[1]
    idx_g2shells[1] = 1

    igl = 1
    for ig = 2:Ng
        if G2_sorted[ig] > (G2_sorted[ig-1] + eps8)
            igl = igl + 1
            G2_shells[igl] = G2_sorted[ig]
        end
        idx_g2shells[ig] = igl
    end

    return G2_shells, idx_g2shells

end

function _flip_fft( mm, S )
    if mm > S/2
        return mm - S
    else
        return mm
    end
end

import Base: show
"""
Display some information about `gvec::GVectors`.
"""
function show( io::IO, gvec::GVectors )
    Ng = gvec.Ng
    G = gvec.G
    G2 = gvec.G2
    
    @printf(io, "\n")
    @printf(io, "                                    --------\n")
    @printf(io, "                                    GVectors\n")
    @printf(io, "                                    --------\n")
    @printf(io, "\n")
    @printf(io, "Ng = %12d\n", Ng)
    @printf(io, "\n")
    for ig = 1:3
        @printf(io, "%8d [%18.10f,%18.10f,%18.10f] : %18.10f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])        
    end
    @printf(io, " ....... \n")
    for ig = Ng-3:Ng
        @printf(io, "%8d [%18.10f.%18.10f,%18.10f] : %18.10f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])
    end
    @printf(io, "\n")
    @printf(io, "Max G2 = %18.10f\n", maximum(G2))
end
show( gvec::GVectors ) = show( stdout, gvec )
