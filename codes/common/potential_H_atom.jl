function pot_H_atom( atoms::Atoms, grid )
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

function pot_Hps_HGH( atoms::Atoms, grid )
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