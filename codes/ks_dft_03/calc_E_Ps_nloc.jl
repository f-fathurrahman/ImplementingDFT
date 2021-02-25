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