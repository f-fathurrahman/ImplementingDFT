function op_V_Ps_nloc( Ham::Hamiltonian, psi::Array{Float64,2} )
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

    for ist = 1:Nstates
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
                    for ip = 1:Npoints
                        Vpsi[ip,ist] = Vpsi[ip,ist] + hij*betaNL[ip,ibeta]*betaNL_psi[ist,jbeta]
                    end
                end # iprj
                end # jprj
            end # m
            end # l
        end
    end
    return Vpsi
end