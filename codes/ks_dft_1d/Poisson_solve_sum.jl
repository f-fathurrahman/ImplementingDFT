# a is in bohr, softening parameter of soft Coulomb potential
function Poisson_solve_sum!( grid,
    rho::Vector{Float64}, V::Vector{Float64}; a = 1.0
)
    xgrid = grid.x
    hx = grid.hx
    Npoints = size(rho,1)
    fill!(V, 0.0)
    for ip in 1:Npoints
        xi = xgrid[ip]
        for jp in 1:Npoints
            xj = xgrid[jp]
            dr = sqrt( (xi - xj)^2 + a^2 )
            V[ip] += rho[jp]/dr
        end
        V[ip] = V[ip]*hx
    end
    return
end 

