using Printf
using LinearAlgebra
using SparseArrays
using IncompleteLU
using AlgebraicMultigrid

include("../3d/FD3dGrid.jl")
include("../3d/build_nabla2_matrix.jl")
include("../supporting_functions.jl")
include("Poisson_solve_CG.jl")
include("Poisson_solve_PCG.jl")

function test_main( NN::Array{Int64} )
    AA = [0.0, 0.0, 0.0]
    BB = [16.0, 16.0, 16.0]

    # Initialize grid
    FD = FD3dGrid( NN, AA, BB )

    # Box dimensions
    Lx = BB[1] - AA[1]
    Ly = BB[2] - AA[2]
    Lz = BB[3] - AA[3]

    # Center of the box
    x0 = Lx/2.0
    y0 = Ly/2.0
    z0 = Lz/2.0

    # Parameters for two gaussian functions
    sigma1 = 0.75
    sigma2 = 0.50
    Npoints = FD.Nx * FD.Ny * FD.Nz

    rho = zeros(Float64, Npoints)
    phi = zeros(Float64, Npoints)

    # Initialization of charge density
    dr = zeros(Float64,3)
    for ip in 1:Npoints
        dr[1] = FD.r[1,ip] - x0
        dr[2] = FD.r[2,ip] - y0
        dr[3] = FD.r[3,ip] - z0
        r = norm(dr)
        rho[ip] = exp( -r^2 / (2.0*sigma2^2) ) / (2.0*pi*sigma2^2)^1.5 -
                  exp( -r^2 / (2.0*sigma1^2) ) / (2.0*pi*sigma1^2)^1.5
    end

    deltaV = FD.hx * FD.hy * FD.hz

    println("Building Laplacian3d")
    Laplacian3d = build_nabla2_matrix( FD, func_1d=build_D2_matrix_9pt )

    println("Building preconditioner")
    #prec = ilu(Laplacian3d)
    #prec = ilu(Laplacian3d, τ = 0.001)
    prec = aspreconditioner(ruge_stuben(Laplacian3d))
    #prec = aspreconditioner(smoothed_aggregation(Laplacian3d))

    @printf("Test norm charge: %18.10f\n", sum(rho)*deltaV)
    print("Solving Poisson equation:\n")

    #phi = Poisson_solve_CG( Laplacian3d, -4*pi*rho, 1000, verbose=true, TOL=1e-10 )
    phi = Poisson_solve_PCG( Laplacian3d, prec, -4*pi*rho, 1000, verbose=true, TOL=1e-10 )

    # Calculation of Hartree energy
    Unum = 0.5*sum( rho .* phi ) * deltaV
    Uana = ( ( 1.0/sigma1 + 1.0/sigma2 ) / 2.0 - sqrt(2.0)/sqrt(sigma1^2 + sigma2^2) ) / sqrt(pi)
    @printf("Numeric  = %18.10f\n", Unum)
    @printf("Uana     = %18.10f\n", Uana)
    @printf("abs diff = %18.10e\n", abs(Unum-Uana))
end

test_main([64,64,64])

