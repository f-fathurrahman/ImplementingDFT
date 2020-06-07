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
    Ncenters = 2
    centers = zeros(3,Ncenters)
    centers[:,1] = [ 0.0, 1.0, 1.0]
    centers[:,2] = [ 0.0, 1.0, 1.0]

    # Parameters for two gaussian functions
    σ = 1.0
    Npoints = grid.Npoints

    rho = zeros(Float64, Npoints)
    phi = zeros(Float64, Npoints)

    phi_analytic = zeros(Float64, Npoints)
    
    # Initialization of charge density
    dr = zeros(Float64,3)
    #nrmfct = (2*pi*σ^2)^1.5
    nrmfct = 1.0
    for ip in 1:Npoints
        for ia in 1:Ncenters
            dr[1] = grid.r[1,ip] - centers[1,ia]
            dr[2] = grid.r[2,ip] - centers[2,ia]
            dr[3] = grid.r[3,ip] - centers[3,ia]
            r = sqrt(dr[1]^2 + dr[2]^2 + dr[3]^2)
            rho[ip] = rho[ip] + exp( -r^2 / (2.0*σ^2) ) / nrmfct
            phi_analytic[ip] = phi_analytic[ip] + (2*pi*σ^2)^1.5 * erf(r/(sqrt(2)*σ))/r / nrmfct
        end
        println(rho[ip])
    end

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
    rho = reshape(rho, (NN[1],NN[2],NN[3]))
    phi_analytic = reshape(phi_analytic, (NN[1],NN[2],NN[3]))

    ix = round(Int64,NN[1])
    iy = round(Int64,NN[2])
    iz = round(Int64,NN[3])
    for ix in 1:NN[1]
        #@printf("%18.10f %18.10f %18.10f\n", grid.x[ix], phi[ix,iy,iz], phi_analytic[ix,iy,iz])
        @printf("%18.10f %18.10f %18.10f\n", grid.x[ix], rho[ix,iy,iz], phi_analytic[ix,iy,iz])
    end
    #for iy in 1:NN[2]
    #    @printf("%18.10f %18.10f %18.10f\n", grid.y[iy], phi[ix,iy,iz], phi_analytic[ix,iy,iz])
    #end
    #for iz in 1:NN[3]
    #    @printf("%18.10f %18.10f %18.10f\n", grid.z[iz], phi[ix,iy,iz], phi_analytic[ix,iy,iz])
    #end

    @printf("Numeric  = %18.10f\n", Unum)
    @printf("Uana     = %18.10f\n", Uana)
    @printf("abs diff = %18.10e\n", abs(Unum-Uana))

    @printf("integ_phi   = %18.10f\n", integ_phi)
    @printf("integ_phi_a = %18.10f\n", integ_phi_a)
    @printf("MAE         = %18.10e\n", abs(integ_phi-integ_phi_a)/Npoints)

     @printf("Test norm charge: %18.10f\n", sum(rho)*dVol)

end

test_main([42,42,42])

