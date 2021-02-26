function calc_energies!( Ham::Hamiltonian, psis::Vector{Array{Float64,2}} )

    dVol = Ham.grid.dVol
    Nspin = size(psis,1)

    # Orbitals contribution
    E_kin = 0.0
    E_Ps_nloc = 0.0
    for i in 1:Nspin
        Ham.ispin = i # not really needed (not accessed in calc_E_kin and calc_E_Ps_nloc)
        E_kin = E_kin + calc_E_kin( Ham, psis[i] )
        if Ham.pspotNL.NbetaNL > 0
            E_Ps_nloc = E_Ps_nloc + calc_E_Ps_nloc( Ham, psis[i] )
        end
    end
    Ham.energies.Kinetic = E_kin
    Ham.energies.Ps_nloc = E_Ps_nloc

    Rhoe_tot = dropdims(sum(Ham.rhoe,dims=2),dims=2)

    Ham.energies.Ps_loc = dot( Ham.V_Ps_loc, Rhoe_tot )*dVol
    
    Ham.energies.Hartree = 0.5*dot( Ham.V_Hartree, Rhoe_tot )*dVol
    
    epsxc = calc_epsxc_VWN( Ham.xc_calc, Ham.rhoe )
    Ham.energies.XC = dot( epsxc, Rhoe_tot ) * dVol

    return
end

function calc_energies(Ham, psi)
    calc_energies!(Ham, psi)
    return Ham.energies
end