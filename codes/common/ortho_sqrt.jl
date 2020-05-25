function ortho_sqrt( psi )
    Udagger = inv(sqrt(psi'*psi))
    return psi*Udagger
end

function ortho_sqrt!( psi )
    Udagger = inv(sqrt(psi'*psi))
    psi[:,:] = psi*Udagger
end

function ortho_sqrt( psi, dVol )
    Udagger = inv(sqrt(psi'*psi*dVol))
    return psi*Udagger
end

function ortho_sqrt!( psi, dVol )
    Udagger = inv(sqrt(psi'*psi*dVol))
    psi[:,:] = psi*Udagger
end