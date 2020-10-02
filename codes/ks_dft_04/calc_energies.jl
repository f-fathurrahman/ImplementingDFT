function calc_E_NN( atoms::Atoms, pspots::Vector{PsPot_GTH} )

    r = atoms.positions
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    Zvals = get_Zvals(pspots) # a safer way

    if atoms.pbc == (true,true,true)
        return calc_E_NN( atoms.LatVecs, atoms, Zvals )
    end

    E_NN = 0.0
    for ia in 1:Natoms
        for ja in ia+1:Natoms
            #
            isp = atm2species[ia]
            jsp = atm2species[ja]
            #
            dx = r[1,ja] - r[1,ia]
            dy = r[2,ja] - r[2,ia]
            dz = r[3,ja] - r[3,ia]
            #
            r_ij = sqrt(dx^2 + dy^2 + dz^2)
            #
            E_NN = E_NN + Zvals[isp]*Zvals[jsp]/r_ij
        end
    end

    return E_NN
end

function calc_E_Ps_nloc( Ham::Hamiltonian, psi::Array{Float64,2} )

    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Natoms = Ham.atoms.Natoms
    atm2species = Ham.atoms.atm2species
    pspots = Ham.pspots
    prj2beta = Ham.pspotNL.prj2beta
    NbetaNL = Ham.pspotNL.NbetaNL

    dVol = Ham.grid.dVol
    betaNL_psi = psi' * Ham.pspotNL.betaNL * dVol

    E_Ps_nloc = 0.0
    for ist = 1:Nstates
        enl1 = 0.0
        for ia = 1:Natoms
            isp = atm2species[ia]
            psp = pspots[isp]
            for l = 0:psp.lmax, m = -l:l
                for iprj = 1:psp.Nproj_l[l+1], jprj = 1:psp.Nproj_l[l+1]
                    ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                    jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                    hij = psp.h[l+1,iprj,jprj]
                    enl1 = enl1 + hij * betaNL_psi[ist,ibeta] * betaNL_psi[ist,jbeta]
                end # jprj, iprj
            end # m, l
        end
        E_Ps_nloc = E_Ps_nloc + Focc[ist]*enl1
    end

    return E_Ps_nloc

end

function calc_E_kin( Ham, psi::Array{Float64,2} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    E_kin = 0.0
    dVol = Ham.grid.dVol
    nabla2psi = zeros(Float64,Nbasis)
    Focc = Ham.electrons.Focc
    for ist in 1:Nstates
        @views nabla2psi = -0.5*Ham.Laplacian*psi[:,ist]
        @views E_kin = E_kin + Focc[ist]*dot( psi[:,ist], nabla2psi[:] )*dVol
    end
    return E_kin
end

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