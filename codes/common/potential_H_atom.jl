function pot_H_atom( atoms::Atoms, grid )
    @assert atoms.pbc == (false,false,false)
    @assert atoms.Nspecies == 1
    Npoints = grid.Npoints
    Vpot = zeros(Npoints)
    Natoms = atoms.Natoms
    atpos = atoms.positions
    for ia in 1:Natoms, ip in 1:Npoints
        dx = grid.r[1,ip] - atpos[1,ia]
        dy = grid.r[2,ip] - atpos[2,ia]
        dz = grid.r[3,ip] - atpos[3,ia]
        Vpot[ip] = Vpot[ip] - 1.0/sqrt(dx^2 + dy^2 + dz^2)
    end
    return Vpot
end

function pot_H_atom_G( atoms::Atoms, grid, gvec::GVectors )
    @assert atoms.pbc == (true,true,true)
    @assert atoms.Nspecies == 1
    Npoints = grid.Npoints
    CellVolume = grid.Lx * grid.Ly * grid.Lz
    Nx = grid.Nx; Ny = grid.Ny; Nz = grid.Nz
    G2 = gvec.G2
    Ng = length(G2)
    if typeof(grid) == LF3dGrid
        # periodic LF
        shifts = [0.5*grid.hx, 0.5*grid.hy, 0.5*grid.hz]
        strf = calc_strfact_shifted( atoms, gvec, shifts )
    else
        # periodic FD    
        strf = calc_strfact( atoms, gvec )
    end
    Vpot = zeros(Float64, Npoints)
    Vg = zeros(ComplexF64, Nx,Ny,Nz)
    isp = 1
    Zval = 1.0
    prefactor = -4*pi*Zval
    for ig = 2:Ng
        Vg[ig] = strf[ig,isp] * prefactor/G2[ig]
    end
    #
    ifft!(Vg)
    @views Vpot[:] = Vpot[:] + real( Vg[:] ) * Npoints / CellVolume
    return Vpot
end



function pot_Hps_GTH( atoms::Atoms, grid )
    @assert atoms.pbc == (false,false,false)
    @assert atoms.Nspecies == 1    
    Npoints = grid.Npoints
    Vpot = zeros( Float64, Npoints )
    # Parameters
    Zval = 1
    rloc = 0.2
    C1 = -4.0663326
    C2 = 0.6678322

    Natoms = atoms.Natoms
    atpos = atoms.positions

    for ia in 1:Natoms, ip in 1:Npoints
        dx2 = ( grid.r[1,ip] - atpos[1,ia] )^2
        dy2 = ( grid.r[2,ip] - atpos[2,ia] )^2
        dz2 = ( grid.r[3,ip] - atpos[3,ia] )^2
        r = sqrt(dx2 + dy2 + dz2)
        if r < eps()
            Vpot[ip] = Vpot[ip] - 2*Zval/(sqrt(2*pi)*rloc) + C1
        else
            rrloc = r/rloc
            Vpot[ip] = Vpot[ip] - Zval/r * erf( r/(sqrt(2.0)*rloc) ) +
                     (C1 + C2*rrloc^2)*exp(-0.5*(rrloc)^2)
        end
    end
    return Vpot
end

function pot_Hps_GTH_G( atoms::Atoms, grid, gvec::GVectors )
    @assert atoms.pbc == (true,true,true)
    @assert atoms.Nspecies == 1
    Npoints = grid.Npoints
    CellVolume = grid.Lx * grid.Ly * grid.Lz
    Nx = grid.Nx; Ny = grid.Ny; Nz = grid.Nz
    G2 = gvec.G2
    Ng = length(G2)
    if typeof(grid) == LF3dGrid
        # periodic LF
        shifts = [0.5*grid.hx, 0.5*grid.hy, 0.5*grid.hz]
        strf = calc_strfact_shifted( atoms, gvec, shifts )
    else
        # periodic FD    
        strf = calc_strfact( atoms, gvec )
    end
    Vpot = zeros(Float64, Npoints)
    Vg = zeros(ComplexF64, Nx,Ny,Nz)
    isp = 1


    # Parameters
    Zval = 1
    rloc = 0.2
    C1 = -4.0663326
    C2 = 0.6678322
    pre1 = -4*pi*Zval
    pre2 = sqrt(8*pi^3)*rloc^3

    # Limiting case G=(0,0,0)
    Vg[1] = 2*pi*Zval*rloc^2 + (2*pi)^1.5 * rloc^3 * (C1 + 3.0*C2)
    Vg[1] = strf[1,isp] * Vg[1]
    for ig in 2:Ng
        Gr = sqrt(G2[ig])*rloc
        expGr2 = exp(-0.5*Gr^2)
        Vg[ig] = pre1/G2[ig]*expGr2 + pre2*expGr2 * (C1 + C2*(3.0 - Gr^2))
        Vg[ig] = strf[ig,isp] * Vg[ig]
    end
    ifft!(Vg)
    @views Vpot[:] = Vpot[:] + real( Vg[:] ) * Npoints / CellVolume
    return Vpot
end