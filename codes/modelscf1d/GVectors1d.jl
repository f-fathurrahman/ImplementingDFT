struct GVectors1d
    Ng::Int64
    G::Vector{Float64}
    G2::Vector{Float64}
    idx_g2r::Vector{Int64}
end

# Nx: sampling points
# Lx: lattice vector (1d)
function GVectors1d(Lx::Float64, Nx::Int64)
    Ng = Nx
    G  = zeros(Float64, Ng)
    G2 = zeros(Float64, Ng)
    idx_g2r = zeros(Int64, Ng)
    #
    ig = 0
    ip = 0
    for i in 0:Nx-1
        ip = ip + 1
        ig = ig + 1 # unnecessary?
        gi = _flip_fft(i, Nx)
        G[ig] = 2π*gi/Lx
        G2[ig] = G[ig]^2
        idx_g2r[ig] = ip # unnecessary
    end
    return GVectors1d(Ng, G, G2, idx_g2r)
end

function _flip_fft( mm, S )
    if mm > S/2
        return mm - S
    else
        return mm
    end
end