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
