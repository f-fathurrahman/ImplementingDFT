function ortho_sqrt( psi )
    Udagger = inv(sqrt(psi'*psi))
    return psi*Udagger
end

function ortho_sqrt!( psi )
    Udagger = inv(sqrt(psi'*psi))
    @views psi[:,:] = psi[:,:]*Udagger
end

function ortho_sqrt( psi, dVol )
    Udagger = inv(sqrt(psi'*psi*dVol))
    return psi*Udagger
end

function ortho_sqrt!( psi, dVol )
    Udagger = inv(sqrt(psi'*psi*dVol))
    @views psi[:,:] = psi[:,:]*Udagger
end