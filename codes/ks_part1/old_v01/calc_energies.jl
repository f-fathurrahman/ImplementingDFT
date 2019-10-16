mutable struct Energies
    Kinetic::Float64
    Ps_loc::Float64
    Hartree::Float64
    XC::Float64
end

import Base: sum
function sum( ene::Energies )
    return ene.Kinetic + ene.Ps_loc + ene.Hartree + ene.XC
end

function calc_E_kin( Ham, psi::Array{Float64,2} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    E_kin = 0.0
    nabla2psi = zeros(Float64,Nbasis)
    dVol = Ham.fdgrid.dVol
    # Assumption: Focc = 2 for all states
    for ist in 1:Nstates
        @views nabla2psi = -0.5*Ham.Laplacian*psi[:,ist]
        E_kin = E_kin + 2.0*dot( psi[:,ist], nabla2psi[:] )*dVol
    end
    return E_kin
end

function calc_energies( Ham::Hamiltonian, psi::Array{Float64,2} )
    dVol = Ham.fdgrid.dVol
    E_kin = calc_E_kin( Ham, psi )
    E_Ps_loc = sum( Ham.V_Ps_loc .* Ham.rhoe )*dVol
    E_Hartree = 0.5*sum( Ham.V_Hartree .* Ham.rhoe )*dVol
    E_XC = sum( excVWN(Ham.rhoe) .* Ham.rhoe )*dVol

    return Energies(E_kin, E_Ps_loc, E_Hartree, E_XC)
end
