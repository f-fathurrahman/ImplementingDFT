function calc_energies!( Ham::Hamiltonian, psi::Array{Float64,2} )

    dVol = Ham.grid.dVol

    Ham.energies.Kinetic = calc_E_kin( Ham, psi )
    Ham.energies.Ps_loc = sum( Ham.V_Ps_loc .* Ham.rhoe )*dVol
    if Ham.pspotNL.NbetaNL > 0
        Ham.energies.Ps_nloc = calc_E_Ps_nloc( Ham, psi )
    end
    Ham.energies.Hartree = 0.5*sum( Ham.V_Hartree .* Ham.rhoe )*dVol
    Ham.energies.XC = sum( excVWN(Ham.rhoe) .* Ham.rhoe )*dVol

    return
end

function calc_energies(Ham, psi)
    calc_energies!(Ham, psi)
    return Ham.energies
end