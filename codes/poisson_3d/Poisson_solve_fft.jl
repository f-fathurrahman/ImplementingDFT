struct PoissonSolverFFT
    pbc::Tuple{Bool,Bool,Bool} # in case FFT is also used for 
    gvec::GVectors
end

# FIXME: Currently only for periodic case
# FIXME: also need to allocate FFT plan?
function PoissonSolverFFT(grid)
    @assert grid.pbc == (true,true,true)
    gvec = GVectors(grid)
    return PoissonSolverFFT( (true,true,true), gvec)
end

function Poisson_solve_fft( grid, gvec::GVectors, rho::Vector{Float64} )
    
    Npoints = grid.Npoints

    Nx = grid.Nx
    Ny = grid.Ny
    Nz = grid.Nz

    ctmp = zeros(ComplexF64,Nx,Ny,Nz)
    ctmp[:] = rho[:]

    # to reciprocal space
    fft!(ctmp)

    ctmp[1] = 0.0 + im*0.0
    for ip in 2:Npoints
        ctmp[ip] = 4.0*pi*ctmp[ip] / gvec.G2[ip]
    end

    # to real space
    ifft!(ctmp)

    return reshape(real(ctmp),Npoints)

end

function Poisson_solve_fft( psolver::PoissonSolverFFT, grid, rho )
    return Poisson_solve_fft( grid, psolver.gvec, rho )
end

function Poisson_solve( psolver::PoissonSolverFFT, grid, rho )
    return Poisson_solve_fft( psolver, grid, rho )
end
