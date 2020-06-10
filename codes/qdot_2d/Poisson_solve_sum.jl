function Poisson_solve_sum( grid::Union{FD2dGrid,LF2dGrid}, rho )
    
    Npoints = grid.Npoints
    V = zeros(Float64,Npoints)
    
    Nx = grid.Nx
    Ny = grid.Nx
    r = grid.r
    dVol = grid.dVol

    ri = zeros(Float64,2)
    rj = zeros(Float64,2)

    for ip in 1:Npoints

        ri[1] = r[1,ip]
        ri[2] = r[2,ip]
        
        for jp in 1:Npoints
            
            rj[1] = r[1,jp]
            rj[2] = r[2,jp]

            if ip == jp
                V[ip] = V[ip] + 2*sqrt(pi)*rho[ip]/sqrt(dVol)
            else
                dr = ( ri[1] - rj[1] )^2 + ( ri[2] - rj[2] )^2
                V[ip] = V[ip] + rho[jp]/sqrt(dr)
            end
        end

        V[ip] = V[ip]*dVol
    
    end

    return V
end 