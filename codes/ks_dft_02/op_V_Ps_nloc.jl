function op_V_Ps_nloc!( Ham::Hamiltonian, psi, Vpsi )
    Nstates = size(psi,2)

    atoms = Ham.atoms
    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    pspots = Ham.pspots
    prj2beta = Ham.pspotNL.prj2beta
    betaNL = Ham.pspotNL.betaNL

    dVol = Ham.grid.dVol
    betaNL_psi = psi' * Ham.pspotNL.betaNL *dVol
    
    Npoints = Ham.grid.Npoints

    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        for l = 0:psp.lmax, m = -l:l
            for iprj = 1:psp.Nproj_l[l+1]
                ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                for jprj = 1:psp.Nproj_l[l+1]
                    jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                    hij = psp.h[l+1,iprj,jprj]
                    for ist in 1:Nstates
                        cc = betaNL_psi[ist,jbeta]*hij
                        # should only loop over limited points (for which the projectors are nonzero)
                        for ip in 1:Npoints
                            Vpsi[ip,ist] = Vpsi[ip,ist] + betaNL[ip,ibeta]*cc
                        end
                    end
                end # iprj
            end # jprj
        end # m, l
    end
    return
end


function op_V_Ps_nloc( Ham::Hamiltonian, psi )
    Nstates = size(psi,2)

    atoms = Ham.atoms
    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    pspots = Ham.pspots
    prj2beta = Ham.pspotNL.prj2beta
    betaNL = Ham.pspotNL.betaNL

    dVol = Ham.grid.dVol
    betaNL_psi = psi' * Ham.pspotNL.betaNL *dVol
    
    Npoints = Ham.grid.Npoints
    Vpsi = zeros(Float64,Npoints,Nstates)

    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        for l = 0:psp.lmax, m = -l:l
            for iprj = 1:psp.Nproj_l[l+1]
                ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                for jprj = 1:psp.Nproj_l[l+1]
                    jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                    hij = psp.h[l+1,iprj,jprj]
                    for ist in 1:Nstates
                        cc = betaNL_psi[ist,jbeta]*hij
                        for ip in 1:Npoints
                            Vpsi[ip,ist] = Vpsi[ip,ist] + betaNL[ip,ibeta]*cc
                        end
                    end
                end # iprj
            end # jprj
        end # m, l
    end
    return Vpsi
end