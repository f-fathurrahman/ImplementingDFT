function op_V_Ps_nloc( Ham::Hamiltonian, psi )
    Nstates = size(psi,2)

    atoms = Ham.atoms
    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    pspots = Ham.pspots
    prj2beta = Ham.pspotNL.prj2beta
    betaNL = Ham.pspotNL.betaNL
    NbetaNL = Ham.pspotNL.NbetaNL

    dVol = Ham.grid.dVol
    betaNL_psi = zeros(Float64,Nstates,NbetaNL)
    for ibeta in 1:NbetaNL
        betaNL_psi[:,ibeta] .= psi'*betaNL[ibeta]*dVol
    end
    #println("betaNL_psi = ")
    #display(betaNL_psi); println()
    
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
                        for (inz, ip) in enumerate(betaNL[ibeta].nzind)
                            Vpsi[ip,ist] = Vpsi[ip,ist] + betaNL[ibeta].nzval[inz]*cc
                        end
                    end
                end # iprj
            end # jprj
        end # m, l
    end
    return Vpsi
end