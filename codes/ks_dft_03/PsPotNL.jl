struct PsPotNL
    NbetaNL::Int64
    prj2beta::Array{Int64,4}
    betaNL::Array{Float64,2}
end

function PsPotNL()
    # return dummy PsPotNL
    betaNL = zeros(Float64,1,1)
    return PsPotNL( 0, zeros(Int64,1,1,1,1), betaNL )
end

function setup_betaNL!( atoms, grid, pspots, betaNL )
    
    Natoms = atoms.Natoms
    Npoints = grid.Npoints
    atm2species = atoms.atm2species

    ibeta = 0
    dr = zeros(3)
    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        for l = 0:psp.lmax
        for m = -l:l
        for iprj = 1:psp.Nproj_l[l+1]
            ibeta = ibeta + 1
            for ip in 1:Npoints
                dr[1] = grid.r[1,ip] - atoms.positions[1,ia]
                dr[2] = grid.r[2,ip] - atoms.positions[2,ia]
                dr[3] = grid.r[3,ip] - atoms.positions[3,ia]
                drm = sqrt( dr[1]^2 + dr[2]^2 + dr[3]^2 )
                betaNL[ip,ibeta] = Ylm_real(l, m, dr)*eval_proj_R(psp, l, iprj, drm)
            end
        end
        end
        end
    end
    #@printf("\nsizeof betaNL  = %18.10f MiB\n", Base.summarysize(betaNL)/1024/1024)
    return
end

function check_betaNL_norm( grid, pspotNL::PsPotNL )
    betaNL = pspotNL.betaNL
    NbetaNL = size(betaNL,2)
    dVol = grid.dVol
    @printf("\nTest normalization of betaNL (Nx=%d, Ny=%d, Nz=%d)\n", grid.Nx, grid.Ny, grid.Nz)
    for ibeta in 1:NbetaNL
        @views ss = dot(betaNL[:,ibeta], betaNL[:,ibeta])*dVol
        @printf("%5d %18.10f\n", ibeta, ss)
    end
    return
end

function setup_betaNL_periodic!( atoms, grid, pspots, betaNL )
    
    Natoms = atoms.Natoms
    Npoints = grid.Npoints
    atm2species = atoms.atm2species
    LL = [grid.Lx, grid.Ly, grid.Lz]
    ibeta = 0
    dr = zeros(3)
    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        for l = 0:psp.lmax
        for m = -l:l
        for iprj = 1:psp.Nproj_l[l+1]
            ibeta = ibeta + 1
            for ip in 1:Npoints
                @views calc_dr_periodic!( LL, grid.r[:,ip], atoms.positions[:,ia], dr )
                drm = sqrt( dr[1]^2 + dr[2]^2 + dr[3]^2 )
                betaNL[ip,ibeta] = Ylm_real(l, m, dr)*eval_proj_R(psp, l, iprj, drm)
            end
        end
        end
        end
    end

    return
end




function PsPotNL( atoms::Atoms, pspots::Array{PsPot_GTH,1}, grid )

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    atpos = atoms.positions

    # 4: indexed from 0:3
    # 0, 1, 2, 3  -> l indexed
    # 1, 2, 3, 4  -> l + 1

    # -3:3
    # -3, -2, -1, 0, 1, 2, 3  -> m
    #  1,  2,  3, 4, 5, 6, 7  -> 4 + m, lmax = 3 + 1

    # -2, -1, 0, 1, 2  -> m
    #  1,  2, 3, 4, 5  -> 3 + m, lmax = 2 + 1

    prj2beta = Array{Int64}(undef,3,Natoms,4,7)
    prj2beta[:] .= -1   # set to invalid index

    NbetaNL = 0
    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        for l = 0:psp.lmax
            for m = -l:l
                for iprj = 1:psp.Nproj_l[l+1]
                    NbetaNL = NbetaNL + 1
                    prj2beta[iprj,ia,l+1,m+psp.lmax+1] = NbetaNL
                end
            end
        end
    end

    # No nonlocal components
    if NbetaNL == 0
        # return dummy PsPotNL
        betaNL = zeros(Float64,1,1)
        return PsPotNL( 0, zeros(Int64,1,1,1,1), betaNL )
    end

    Npoints = grid.Npoints
    betaNL = zeros(Float64, Npoints, NbetaNL)

    if grid.pbc == (true,true,true)
        setup_betaNL_periodic!( atoms, grid, pspots, betaNL )
    else
        setup_betaNL!( atoms, grid, pspots, betaNL )
    end

    return PsPotNL( NbetaNL, prj2beta, betaNL )

end

