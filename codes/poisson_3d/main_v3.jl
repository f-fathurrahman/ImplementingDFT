using Printf
using LinearAlgebra
using SparseArrays

include("INC_poisson_3d.jl")
include("Poisson_solve_DAGE.jl")

function test_main( NN::Array{Int64} )
    
    AA = [-8.0, -8.0, -8.0]
    BB = [ 8.0,  8.0,  8.0]
    
    grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )

    psolver = PoissonSolverDAGE(grid)

    # Box dimensions
    Lx = BB[1] - AA[1]
    Ly = BB[2] - AA[2]
    Lz = BB[3] - AA[3]

    # Center of the box
    x0 = AA[1] + Lx/2.0
    y0 = AA[2] + Ly/2.0
    z0 = AA[3] + Lz/2.0

    # Parameters for two gaussian functions
    sigma1 = 0.75
    sigma2 = 0.50
    Npoints = grid.Npoints

    rho = zeros(Float64, Npoints)
    phi = zeros(Float64, Npoints)

    # Initialization of charge density
    dr = zeros(Float64,3)
    for ip in 1:Npoints
        dr[1] = grid.r[1,ip] - x0
        dr[2] = grid.r[2,ip] - y0
        dr[3] = grid.r[3,ip] - z0
        r = sqrt(dr[1]^2 + dr[2]^2 + dr[3]^2)
        rho[ip] = exp( -r^2 / (2.0*sigma2^2) ) / (2.0*pi*sigma2^2)^1.5 -
                  exp( -r^2 / (2.0*sigma1^2) ) / (2.0*pi*sigma1^2)^1.5
    end

    dVol = grid.dVol

    @printf("Test norm charge: %18.10f\n", sum(rho)*dVol)
    print("Solving Poisson equation:\n")

    phi = Poisson_solve_DAGE(psolver, grid, rho)

    # Calculation of Hartree energy
    Unum = 0.5*sum( rho .* phi ) * dVol
    Uana = ( ( 1.0/sigma1 + 1.0/sigma2 ) / 2.0 - sqrt(2.0)/sqrt(sigma1^2 + sigma2^2) ) / sqrt(pi)
    @printf("Numeric  = %18.10f\n", Unum)
    @printf("Uana     = %18.10f\n", Uana)
    @printf("abs diff = %18.10e\n", abs(Unum-Uana))

end

test_main([64,64,64])

