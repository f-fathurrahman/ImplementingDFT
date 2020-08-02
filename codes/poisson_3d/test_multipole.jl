using Printf
using LinearAlgebra
using SparseArrays
using SpecialFunctions: erf

include("INC_poisson_3d.jl")

function test_main( NN::Array{Int64} )
    AA = [-7.0, -7.0, -7.0]
    BB = [ 7.0,  7.0,  7.0]
    
    grid = FD3dGrid( NN, AA, BB )

    # Box dimensions
    Lx = BB[1] - AA[1]
    Ly = BB[2] - AA[2]
    Lz = BB[3] - AA[3]

    # Center of the box
    x0 = AA[1] + Lx/2.0
    y0 = AA[2] + Ly/2.0
    z0 = AA[3] + Lz/2.0

    println("x0 = ", x0)
    println("y0 = ", y0)
    println("z0 = ", z0)

    # Parameters for two gaussian functions
    σ = 1.0
    Npoints = grid.Npoints

    rho = zeros(Float64, Npoints)
    phi = zeros(Float64, Npoints)

    phi_analytic = zeros(Float64, Npoints)
    # Initialization of charge density
    nrmfct = (2*pi*σ^2)^1.5
    dr = zeros(Float64,3)
    SMALL = eps()
    for ip in 1:Npoints
        dr[1] = grid.r[1,ip] - x0 + SMALL
        dr[2] = grid.r[2,ip] - y0 + SMALL
        dr[3] = grid.r[3,ip] - z0 + SMALL
        r = sqrt(dr[1]^2 + dr[2]^2 + dr[3]^2)
        rho[ip] = exp( -r^2 / (2.0*σ^2) )/nrmfct
        phi_analytic[ip] = (2*pi*σ^2)^1.5 * erf(r/(sqrt(2)*σ))/r/nrmfct
    end

    dVol = grid.dVol

    println("sum phi_analytic = ", sum(phi_analytic))
    println("integ rho = ", sum(rho)*dVol)

    lmax = 4
    lmmax = (lmax+1)^2
    Q_lm = zeros(Float64,lmmax)

    calc_multipole_moment!( grid, rho, Q_lm )

    println("Multipole moment: ")
    for i in 1:lmmax
        @printf("%3d %18.10f\n", i, Q_lm[i])
    end

    V_x0 = zeros(Float64, grid.Ny, grid.Nz)
    V_xN = zeros(Float64, grid.Ny, grid.Nz)

    V_y0 = zeros(Float64, grid.Nx, grid.Nz)
    V_yN = zeros(Float64, grid.Nx, grid.Nz)

    V_z0 = zeros(Float64, grid.Nx, grid.Ny)
    V_zN = zeros(Float64, grid.Nx, grid.Ny)

    @time set_bc_isolated!( grid, Q_lm, V_x0, V_xN, V_y0, V_yN, V_z0, V_zN )

end

test_main([50,50,50])