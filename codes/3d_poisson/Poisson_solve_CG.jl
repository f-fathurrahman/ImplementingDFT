function Poisson_solve_CG( Lmat::SparseMatrixCSC{Float64,Int64},
                           rho::Array{Float64,1}, NiterMax::Int64;
                           verbose=false, TOL=5.e-10 )
  #
  Npoints = size(rho)[1]
  phi = zeros( Float64, Npoints ) # XXX or use some starting guess
  #
  r = zeros( Float64, Npoints )
  p = zeros( Float64, Npoints )
  #
  nabla2_phi = Lmat*phi
  for ip = 1 : Npoints
    r[ip] = rho[ip] - nabla2_phi[ip]
    p[ip] = r[ip]
  end
  rsold = dot( r, r )

  for iter = 1 : NiterMax
    #
    nabla2_phi = Lmat*p
    #
    alpha = rsold/dot( p, nabla2_phi )
    #
    phi = phi + alpha * p
    r = r - alpha * nabla2_phi
    #
    rsnew = dot( r, r )
    deltars = rsold - rsnew
    if verbose
      @printf("%8d %20.10f\n", iter, sqrt(rsnew))
    end
    #
    if sqrt(rsnew) < TOL
      if verbose
        @printf("Convergence achieved in Poisson_solve_CG: N, iter: %d %d\n", Npoints, iter)
      end
      break
    end
    #
    p = r + (rsnew/rsold) * p
    #
    rsold = rsnew
  end
  #
  return phi
  #
end # of function
