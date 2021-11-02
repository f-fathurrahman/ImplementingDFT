# Hamiltonian operator
function op_H( Ham::Hamiltonian, psi )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    Hpsi = zeros(Float64,Nbasis,Nstates)
    #
    #Hpsi = -0.5*Ham.∇2 * psi
    mul!(Hpsi, Ham.∇2, psi)
    lmul!(-0.5, Hpsi)
    #
    if Ham.pspotNL.NbetaNL > 0
        Vnlpsi = op_V_Ps_nloc(Ham, psi)
        for ist in 1:Nstates, ip in 1:Nbasis
            Hpsi[ip,ist] = Hpsi[ip,ist] + ( Ham.V_Ps_loc[ip] +
                Ham.V_Hartree[ip] + Ham.V_XC[ip] ) * psi[ip,ist] + Vnlpsi[ip,ist]
        end
    else
        for ist in 1:Nstates, ip in 1:Nbasis
            Hpsi[ip,ist] = Hpsi[ip,ist] + ( Ham.V_Ps_loc[ip] +
                Ham.V_Hartree[ip] + Ham.V_XC[ip] ) * psi[ip,ist]
        end
    end
    return Hpsi
end

import Base: *
function *( Ham::Hamiltonian, psi )
    return op_H(Ham, psi)
end