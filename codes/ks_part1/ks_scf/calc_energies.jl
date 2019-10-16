function calc_E_kin( Ham, psi::Array{Float64,2} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    E_kin = 0.0
    dVol = Ham.fdgrid.dVol
    nabla2psi = zeros(Float64,Nbasis)
    for ist in 1:Nstates
        @views nabla2psi = -0.5*Ham.Laplacian*psi[:,ist]
        E_kin = E_kin + 2.0*dot( psi[:,ist], nabla2psi[:] )*dVol
    end
    return E_kin
end

function calc_energies!( Ham::Hamiltonian, psi::Array{Float64,2} )
    
    dVol = Ham.fdgrid.dVol

    Ham.energies.Kinetic = calc_E_kin( Ham, psi )
    Ham.energies.Ps_loc = sum( Ham.V_Ps_loc .* Ham.rhoe )*dVol
    Ham.energies.Hartree = 0.5*sum( Ham.V_Hartree .* Ham.rhoe )*dVol
    Ham.energies.XC = sum( excVWN(Ham.rhoe) .* Ham.rhoe )*dVol

    return !
end