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

    # calculate E_NL
    E_Ps_nloc = 0.0

    dVol = Ham.grid.dVol
    betaNL_psi = psi' * Ham.pspotNL.betaNL * dVol

    for ist = 1:Nstates
        enl1 = 0.0
        for ia = 1:Natoms
            isp = atm2species[ia]
            psp = pspots[isp]
            for l = 0:psp.lmax
            for m = -l:l
            for iprj = 1:psp.Nproj_l[l+1]
            for jprj = 1:psp.Nproj_l[l+1]
                ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                hij = psp.h[l+1,iprj,jprj]
                enl1 = enl1 + hij * betaNL_psi[ist,ibeta] * betaNL_psi[ist,jbeta]
            end
            end
            end # m
            end # l
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
        E_kin = E_kin + Focc[ist]*dot( psi[:,ist], nabla2psi[:] )*dVol
    end
    return E_kin
end

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
