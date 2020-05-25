function Poisson_solve_PCG( Lmat::SparseMatrixCSC{Float64,Int64},
                            prec,
                            rho::Array{Float64,1}, NiterMax::Int64;
                            verbose=false, TOL=5.e-10 )
    #
    Npoints = size(rho,1)
    phi = zeros( Float64, Npoints ) # XXX or use some starting guess
    #
    r = zeros( Float64, Npoints )
    p = zeros( Float64, Npoints )
    z = zeros( Float64, Npoints )
    #
    nabla2_phi = Lmat*phi
    r = rho - nabla2_phi
    z = copy(r)
    ldiv!(prec, z)
    p = copy(z)

    rsold = dot( r, z )

    for iter = 1 : NiterMax
        #
        nabla2_phi = Lmat*p
        #
        alpha = rsold/dot( p, nabla2_phi )
        #
        phi = phi + alpha * p
        r = r - alpha * nabla2_phi
        z = copy(r)
        ldiv!(prec, z)
        #
        rsnew = dot(z, r)
        deltars = rsold - rsnew

        if verbose
            @printf("%8d %18.10e\n", iter, sqrt(abs(rsnew)))
        end
        #
        if sqrt(abs(rsnew)) < TOL
            if verbose
                @printf("Convergence achieved in Poisson_solve_PCG in iter: %d\n", iter)
            end
            break
        end
        #
        p = z + (rsnew/rsold) * p
        #
        rsold = rsnew
    end
    #
    return phi
    #
end

