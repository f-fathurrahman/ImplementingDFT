function calc_E_kin( Ham, psi::Array{Float64,2} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    E_kin = 0.0
    dVol = Ham.grid.dVol
    nabla2psi = zeros(Float64,Nbasis)
    Focc = Ham.electrons.Focc
    for ist in 1:Nstates
        @views nabla2psi = -0.5*Ham.Laplacian*psi[:,ist]
        E_kin = E_kin + Focc[ist]*dot( psi[:,ist], nabla2psi[:] )*dVol
    end
    return E_kin
end

function calc_energies!( Ham::Hamiltonian, psi::Array{Float64,2} )

    dVol = Ham.grid.dVol

    Ham.energies.Kinetic = calc_E_kin( Ham, psi )
    Ham.energies.Ps_loc = dot( Ham.V_Ps_loc, Ham.rhoe )*dVol
    Ham.energies.Hartree = 0.5*dot( Ham.V_Hartree, Ham.rhoe )*dVol
    Ham.energies.XC = calc_E_xc_2d( Ham.rhoe, dVol=dVol )

    return
end
