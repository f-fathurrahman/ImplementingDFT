function op_H( Ham::Hamiltonian, psi::Matrix{Float64} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    Hpsi = zeros(Float64,Nbasis,Nstates)
    #
    Hpsi = -0.5*Ham.Laplacian * psi
    #
    ispin = Ham.ispin # !!! Need to get active spin index
    #
    if Ham.pspotNL.NbetaNL > 0
        Vnlpsi = op_V_Ps_nloc(Ham, psi)
        for ist in 1:Nstates, ip in 1:Nbasis
            Hpsi[ip,ist] = Hpsi[ip,ist] + ( Ham.V_Ps_loc[ip] +
                Ham.V_Hartree[ip] + Ham.V_XC[ip,ispin] ) * psi[ip,ist] + Vnlpsi[ip,ist]
        end
    else
        for ist in 1:Nstates, ip in 1:Nbasis
            Hpsi[ip,ist] = Hpsi[ip,ist] + ( Ham.V_Ps_loc[ip] +
                Ham.V_Hartree[ip] + Ham.V_XC[ip,ispin] ) * psi[ip,ist]
        end
    end
    return Hpsi
end

import Base: *
function *(Ham, psi)
    return op_H(Ham, psi)
end

function op_H( Ham::Hamiltonian, psis::Vector{Matrix{Float64}} )
    Nspin = size(psis,1)
    Hpsis = Vector{Matrix{Float64}}(undef,Nspin)
    for i in 1:Nspin
        Ham.ispin = i
        Hpsis[i] = op_H( Ham, psis[i] )
    end
end