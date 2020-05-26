using SpecialFunctions: erf, erfcx

struct PoissonSolverDAGE
    N_t::Int64
    t_i::Float64
    t_l::Float64
    t_f::Float64
    w_t::Array{Float64,1}
    x_t::Array{Float64,1}
    F_xs::Array{Float64,3}
    F_ys::Array{Float64,3}
    F_zs::Array{Float64,3}
end

# Generate Gauss-Legendre N-points quadrature formula
function _init_gauss_legendre( x1::Float64, x2::Float64, N::Int64 )

    SMALL = 3.e-11

    m = round(Int64,(N + 1)/2) # The roots are symmetric in the interval, so
    xm = 0.5*(x2 + x1)         # we only have to find half of them.
    xl = 0.5*(x2 - x1)

    # FIXME: Need this ?
    z1 = 100.0  # some initial number to get the while loop enter the first time
    pp = 0.0    # make pp visible outside for loop

    x = zeros(Float64,N)
    w = zeros(Float64,N)

    for i in 1:m # Loop over the desired roots.
        z = cos( pi*(i-0.25)/(N+0.5) )
        # Starting with the above approximation to the ith root, we enter the main loop of
        # refinement by Newton’s method.
        while true
            p1 = 1.0
            p2 = 0.0
            for j in 1:N      # Loop up the recurrence relation to get the
                p3 = p2       # Legendre polynomial evaluated at z.
                p2 = p1
                p1 = ( (2.0*j - 1.0)*z*p2 - (j - 1.0)*p3 )/j
            end 
            # p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
            # by a standard relation involving also p2, the polynomial of one lower order.
            pp = N*(z*p1-p2)/(z*z-1.0)
            z1 = z
            z  = z1 - p1/pp # Newton’s method.
            if abs(z-z1) < SMALL
                break
            end
        end
    
        x[i]     = xm - xl*z                  # Scale the root to the desired interval,
        x[N-i+1] = xm + xl*z                  # and put in its symmetric counterpart.
        w[i]     = 2.0*xl/((1.0-z*z)*pp*pp)   # Compute the weight
        w[N-i+1] = w[i]                       # and its symmetric counterpart.
    end

    return x, w

end 


# Sundholm (JCP 132, 024102, 2010).
# Within two-divided regions:
#     ([t_i,t_l], [t_l,t_f]),
# num_points1, num_points2 quadrature points are made, respectively.
function _init_t_sampling( num_points1, num_points2, t_i, t_l, t_f )

    t_values = zeros(Float64,num_points1 + num_points2)
    w_t = zeros(Float64,num_points1 + num_points2)

    # Linear coord region:  [t_i, t_l]
    x_leg, w_leg = _init_gauss_legendre(t_i, t_l, num_points1)

    for j in 1:num_points1
        t_values[j] = x_leg[j]
        w_t[j]      = w_leg[j]
    end

    # Logarithmic coord region: [t_l, t_f]
    x_leg2, w_leg2 = _init_gauss_legendre( log(t_l), log(t_f), num_points2 )

    # Return the log-coord-partitioned points back to linear t-space.
    s_p = 0.0
    w_p = 0.0
    for j in 1:num_points2
        s_p = x_leg2[j]
        w_p = w_leg2[j]
        x_leg2[j] = exp(s_p)
        w_leg2[j] = w_p * exp(s_p)
    end 

    for j in 1:num_points2
        t_values[num_points1+j] = x_leg2[j]
        w_t[num_points1+j]      = w_leg2[j]
    end

    return t_values, w_t

end


function _compute_F( t, x_bar, h )

    f = 0.0
    SMALL = eps()

    if x_bar < SMALL
        f = sqrt(h) * erf( pi/(2.0*h*t) )
    else 
        z = pi/(2.0*h*t) + im*t*x_bar
        w_iz = erfcx( z )
        f = exp( -t*t*x_bar*x_bar )
        f = f - real( exp( -t*t*x_bar*x_bar - z*z )*w_iz ) 
        f = sqrt(h)*f
    end

    return f
end

function _construct_F( grd::Vector{Float64}, N::Int64, h::Float64, t_size, t_values )
    F_values = zeros(Float64,N,N,t_size)
    for i_t in 1:t_size
        for j in 1:N, i in 1:N
            x_bar = abs(grd[i] - grd[j])
            F_values[i,j,i_t] = _compute_F( t_values[i_t], x_bar, h )
        end
    end
    return F_values
end



function PoissonSolverDAGE(
    grid;
    t_i=0.0, t_l=1.0, t_f=100_000.0,
    num_points1=50, num_points2=50
)
  N_t = num_points1 + num_points2

  Nx = grid.Nx
  Ny = grid.Ny
  Nz = grid.Nz

  F_xs = zeros(Float64,Nx,Nx,N_t)
  F_ys = zeros(Float64,Ny,Ny,N_t)
  F_zs = zeros(Float64,Nz,Nz,N_t)

  x_t = zeros(Float64,N_t)
  w_t = zeros(Float64,N_t)

  x_t, w_t = _init_t_sampling( num_points1, num_points2, t_i, t_l, t_f )

  F_xs = _construct_F( grid.x, grid.Nx, grid.hx, N_t, x_t )
  F_ys = _construct_F( grid.y, grid.Ny, grid.hy, N_t, x_t )
  F_zs = _construct_F( grid.z, grid.Nz, grid.hz, N_t, x_t )

  return PoissonSolverDAGE( N_t, t_i, t_l, t_f, w_t, x_t, F_xs, F_ys, F_zs )

end


function Poisson_solve_DAGE( psolver::PoissonSolverDAGE, grid, Rhoe::Vector{Float64} )

    Nx = grid.Nx
    Ny = grid.Ny
    Nz = grid.Nz
    Npoints = grid.Npoints
    dVol = grid.dVol

    density = reshape(Rhoe,Nx,Ny,Nz)
  
    T_g = zeros(Float64,Nx,Ny,Nz)
    T_g2 = zeros(Float64,Nx,Ny,Nz)
    T_b = zeros(Float64,Nx,Ny,Nz)
    T_b2 = zeros(Float64,Nx,Ny,Nz)

    potential = zeros(Float64,Nx,Ny,Nz)

    F_xs = psolver.F_xs
    F_ys = psolver.F_ys
    F_zs = psolver.F_zs
    N_t = psolver.N_t
    w_t = psolver.w_t
    t_f = psolver.t_f

    for i_t in 1:N_t
        for gg in 1:Nz
            @views T_g[:,:,gg] = F_xs[:,:,i_t] * density[:,:,gg]
            @views T_g2[:,:,gg] = T_g[:,:,gg] * F_ys[:,:,i_t]
        end

        # reorder
        for bb = 1:Ny, gg in 1:Nz
            for ix in 1:Nx
                T_b[ix,gg,bb] = T_g2[ix,bb,gg]
            end
        end

        for bb in 1:Ny
            @views T_b2[:,:,bb] = T_b[:,:,bb] * F_zs[:,:,i_t]
        end

        for k in 1:Nz, j in 1:Ny, i in 1:Nx
            potential[i,j,k] = potential[i,j,k] + w_t[i_t]*T_b2[i,j,k]*2.0/sqrt(pi)
        end

    end

    for ip in 1:Npoints
        potential[ip] = ( pi/(t_f*t_f) ) * density[ip] + potential[ip]
    end
    
    return reshape(potential,Npoints)*sqrt(dVol)

end
