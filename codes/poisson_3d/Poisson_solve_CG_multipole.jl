#
# Adapted from Varga - Computational Nanoscience
#

function Ylm( lm::Int64, x::Float64, y::Float64, z::Float64 )
    
    res = 0.0
    
    # Spherical harmonics    
    # It should be multiplied by r^l*sqrt((2*l+1)/4*pi)
    r = x^2 + y^2 + z^2
    
    if lm == 1
        return 1.0                       # lm=1  (0  0)
    #
    elseif lm == 2
        return y                         # lm=2  (1 -1)
    elseif lm == 3
        return z                         # lm=3  (1  0)
    elseif lm == 4
        return x                         # lm=4  (1  1)
    #
    elseif lm == 5
        return sqrt(3.0)*x*y             # lm=5  (2 -2)
    elseif lm == 6
        return sqrt(3.0)*y*z             # lm=6  (2 -1)
    elseif lm == 7
        return (2*z*z-x*x-y*y)/2.0       # lm=7  (2  0)
    elseif lm == 8
        return sqrt(3.0)*x*z             # lm=8  (2  1)
    elseif lm == 9
        return sqrt(3.0/4.0)*(x*x-y*y)   # lm=9  (2  2)
    #
    elseif lm == 10
        return sqrt(5.0/8.0)*y*(3*x*x-y*y)            # lm=10 (3 -3)
    elseif lm == 11
        return sqrt(15.0)*x*y*z                       # lm=11 (3 -2)    
    elseif lm == 12
        return sqrt(3.0/8.0)*y*(4*z*z-x*x-y*y)        # lm=12 (3 -1)    
    elseif lm == 13
        return z*(2*z*z-3*x*x-3*y*y)/2.0              # lm=13 (3  0)
    elseif lm == 14
        return sqrt(3.0/8.0)*x*(4*z*z-x*x-y*y)        # lm=14 (3  1)
    elseif lm == 15
        return sqrt(15.0/4.0)*z*(x*x-y*y)             # lm=15 (3  2)
    elseif lm == 16
        return sqrt(5.0/8.0)*x*(x*x-3*y*y)            # lm=16 (3  3)
    #
    elseif lm == 17
        return sqrt(35.0)/2.0*x*y*(x^2-y^2)         # lm=17 (4 -4)
    elseif lm == 18
        return sqrt(35.0/8.0)*y*z*(3*x^2-y^2)       # lm=18 (4 -3)
    elseif lm == 19
        return sqrt(5.0)/2.0*x*y*(7*z^2-r^2)        # lm=19 (4 -2)
    elseif lm == 20
        return sqrt(5.0/8.0)*y*(7*z^3-3*z*r^2)      # lm=20 (4 -1)
    elseif lm == 21
        return (35*z^4-30*z^2*r^2+3.0*r^2)/8.0      # lm=21 (4  0)
    elseif lm == 22
        return sqrt(5.0/8.0)*x*(7*z^3-3*z*r^2)      # lm=22 (4  1)
    elseif lm == 23
        return sqrt(5.0)/4.0*(7*z^2-r^2)*(x^2-y^2)  # lm=23 (4  2)
    elseif lm == 24
        return sqrt(35.0/8.0)*z*x*(x^2-3*y^2)       # lm=24 (4  3)
    elseif lm == 25
        return sqrt(35.0)/8.0*(x^4+y^4-6*x^2*y^2)   # lm=25 (4  4)
    #
    else
        return 0.0
    end

end


# Calculate Q_lm
function calc_multipole_moment!( grid, rho, Q_lm ; L_max=4 )
    SMALL = eps()
    Q_lm .= 0.0
    Npoints = grid.Npoints
    dVol = grid.dVol
    for L in 0:L_max
        for lm in (L^2+1):((L+1)^2)
            for ip in 1:Npoints
                x  = grid.r[1,ip]
                y  = grid.r[2,ip]
                z  = grid.r[3,ip]
                r  = sqrt(x*x + y*y + z*z) + SMALL
                xx = x/r
                yy = y/r
                zz = z/r
                Q_lm[lm] = Q_lm[lm] + r^L * Ylm(lm, xx, yy, zz) * rho[ip] * dVol
            end
        end
    end
    return
end

function calc_multipole_potential( Q_lm, L_max, x, y, z )
    SMALL = eps()
    r = sqrt(x*x + y*y + z*z) + SMALL
    xx = x/r
    yy = y/r
    zz = z/r
    V_lm = 0.0
    for L in 0:L_max
        for lm in (L^2+1):((L+1)^2)
            V_lm = V_lm + Ylm(lm, xx, yy, zz)/r^(L+1) * Q_lm[lm]
        end
    end
    return V_lm
end

# Set BC for isolated system
function set_bc_isolated!( grid, Q_lm, V_x0, V_xN, V_y0, V_yN, V_z0, V_zN  )

    idx_xyz2ip = grid.idx_xyz2ip

    # Boundary condition in the x direction
    for k in 1:grid.Nz
        for j in 1:grid.Ny
            #
            y = grid.r[2,idx_xyz2ip[1,j,k]]
            z = grid.r[3,idx_xyz2ip[1,j,k]]
            #
            x = grid.r[1,idx_xyz2ip[1,j,k]] - grid.hx
            V_x0[j,k] = calc_multipole_potential( Q_lm, 4, x, y, z )
            #
            x = grid.r[1,idx_xyz2ip[grid.Nx,j,k]] + grid.hx
            V_xN[j,k] = calc_multipole_potential( Q_lm, 4, x, y, z )
        end
    end

    # Boundary condition in the y direction
    for k in 1:grid.Nz
        for i in 1:grid.Nx
            #
            x = grid.r[1,idx_xyz2ip[i,1,k]]
            z = grid.r[3,idx_xyz2ip[i,1,k]]
            #
            y = grid.r[2,idx_xyz2ip[i,1,k]] - grid.hy
            V_y0[i,k] = calc_multipole_potential( Q_lm, 4, x, y, z )
            #
            y = grid.r[2,idx_xyz2ip[i,grid.Ny,k]] + grid.hy
            V_yN[i,k] = calc_multipole_potential( Q_lm, 4, x, y, z )
        end
    end

    # Boundary condition in the z direction
    for j in 1:grid.Ny
        for i in 1:grid.Nx
            #
            x = grid.r[1,idx_xyz2ip[i,j,1]]
            y = grid.r[2,idx_xyz2ip[i,j,1]]
            #
            z = grid.r[3,idx_xyz2ip[i,j,1]] - grid.hz
            V_z0[i,j] = calc_multipole_potential( Q_lm, 4, x, y, z )
            #
            z = grid.r[3,idx_xyz2ip[i,j,grid.Nz]] + grid.hz
            V_zN[i,j] = calc_multipole_potential( Q_lm, 4, x, y, z )
        end
    end

    return
end