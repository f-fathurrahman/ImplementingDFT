function pot_H_atom( grid; r0=(0.0, 0.0, 0.0) )
    
    Npoints = grid.Npoints
    @assert iseven(Npoints) # Caution: make sure to avoid the singularity

    Vpot = zeros(Npoints)
    for i in 1:Npoints
        dx = grid.r[1,i] - r0[1]
        dy = grid.r[2,i] - r0[2]
        dz = grid.r[3,i] - r0[3]
        Vpot[i] = -1.0/sqrt(dx^2 + dy^2 + dz^2)
    end
    return Vpot
end

function pot_Hps_GTH( grid; r0=(0.0, 0.0, 0.0) )
    Npoints = grid.Npoints
    Vpot = zeros( Float64, Npoints )

    # Parameters
    Zval = 1
    rloc = 0.2
    C1 = -4.0663326
    C2 = 0.6678322

    # TODO Add journal reference
    for ip = 1:Npoints
        dx2 = ( grid.r[1,ip] - r0[1] )^2
        dy2 = ( grid.r[2,ip] - r0[2] )^2
        dz2 = ( grid.r[3,ip] - r0[3] )^2
        r = sqrt(dx2 + dy2 + dz2)
        if r < eps()
            Vpot[ip] = -2*Zval/(sqrt(2*pi)*rloc) + C1
        else
            rrloc = r/rloc
            Vpot[ip] = -Zval/r * erf( r/(sqrt(2.0)*rloc) ) +
                     (C1 + C2*rrloc^2)*exp(-0.5*(rrloc)^2)
        end
    end
    return Vpot
end