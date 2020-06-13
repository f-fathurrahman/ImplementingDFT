# Calculate structure factor
function calc_strfact(
    atpos::Array{Float64,2}, Nspecies::Int64,
    atm2species::Array{Int64,1}, G::Array{Float64,2} 
)
    Ng = size(G)[2]
    Na = size(atpos)[2]
    Sf = zeros(ComplexF64,Ng,Nspecies)
    for ia = 1:Na
        isp = atm2species[ia]
        for ig = 1:Ng
            GX = atpos[1,ia]*G[1,ig] + atpos[2,ia]*G[2,ig] + atpos[3,ia]*G[3,ig]
            Sf[ig,isp] = Sf[ig,isp] + cos(GX) - im*sin(GX)
        end
    end
    return Sf
end

function calc_strfact( atoms::Atoms, gvec::GVectors )
    return calc_strfact( atoms.positions, atoms.Nspecies, atoms.atm2species, gvec.G )
end

#
# Shifted version (for LF3d periodic)
#
function calc_strfact_shifted(
    atpos::Array{Float64,2}, Nspecies::Int64,
    atm2species::Array{Int64,1}, G::Array{Float64,2}, Δr
)
    Ng = size(G)[2]
    Na = size(atpos)[2]
    Sf = zeros(ComplexF64,Ng,Nspecies)
    for ia = 1:Na
        isp = atm2species[ia]
        for ig = 1:Ng
            GX = (atpos[1,ia] - Δr[1]) * G[1,ig] +
                 (atpos[2,ia] - Δr[2]) * G[2,ig] +
                 (atpos[3,ia] - Δr[3]) * G[3,ig]
            Sf[ig,isp] = Sf[ig,isp] + cos(GX) - im*sin(GX)
        end
    end
    return Sf
end

function calc_strfact_shifted( atoms::Atoms, gvec::GVectors, Δr )
    return calc_strfact( atoms.positions, atoms.Nspecies, atoms.atm2species, gvec.G, Δr )
end
