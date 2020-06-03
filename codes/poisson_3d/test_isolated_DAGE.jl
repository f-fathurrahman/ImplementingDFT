using Printf
using LinearAlgebra
using SparseArrays
using AlgebraicMultigrid

using SpecialFunctions: erf

include("INC_poisson_3d.jl")

function test_main( NN::Array{Int64} )
    AA = [-7.0, -7.0, -7.0]
    BB = [ 7.0,  7.0,  7.0]

    #grid = FD3dGrid( NN, AA, BB )
    grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )

    println(grid)

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
    σ = 1.0
    Npoints = grid.Npoints

    rho = zeros(Float64, Npoints)
    phi = zeros(Float64, Npoints)

    phi_analytic = zeros(Float64, Npoints)
    
    # Initialization of charge density
    dr = zeros(Float64,3)
    nrmfct = (2*pi*σ^2)^1.5
    #nrmfct = 1.0
    for ip in 1:Npoints
        dr[1] = grid.r[1,ip] - x0
        dr[2] = grid.r[2,ip] - y0
        dr[3] = grid.r[3,ip] - z0
        r = sqrt(dr[1]^2 + dr[2]^2 + dr[3]^2)
        rho[ip] = exp( -r^2 / (2.0*σ^2) ) / nrmfct
        phi_analytic[ip] = (2*pi*σ^2)^1.5 * erf(r/(sqrt(2)*σ))/r / nrmfct
    end

    println("sum phi_analytic = ", sum(phi_analytic))

    dVol = grid.dVol

    @printf("Test norm charge: %18.10f\n", sum(rho)*dVol)
    print("Solving Poisson equation:\n")

    phi = Poisson_solve_DAGE(psolver, grid, rho)
    println("sum phi = ", sum(phi))

    # Calculation of Hartree energy
    Unum = 0.5*sum( rho .* phi ) * dVol
    Uana = 0.5*sum( rho .* phi_analytic ) * dVol

    integ_phi   = sum( phi ) * dVol
    integ_phi_a = sum( phi_analytic ) * dVol

    phi = reshape(phi, (NN[1],NN[2],NN[3]))
    phi_analytic = reshape(phi_analytic, (NN[1],NN[2],NN[3]))

    ix = NN[1]
    iz = NN[3]
    for iy in 1:NN[2]
        @printf("%18.10f %18.10f %18.10f\n", grid.y[iy], phi[ix,iy,iz], phi_analytic[ix,iy,iz])
    end

    @printf("Numeric  = %18.10f\n", Unum)
    @printf("Uana     = %18.10f\n", Uana)
    @printf("abs diff = %18.10e\n", abs(Unum-Uana))

    @printf("integ_phi   = %18.10f\n", integ_phi)
    @printf("integ_phi_a = %18.10f\n", integ_phi_a)
    @printf("MAE         = %18.10e\n", abs(integ_phi-integ_phi_a)/Npoints)

end

test_main([50,50,50])

