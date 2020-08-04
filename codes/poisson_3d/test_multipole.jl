using Printf
using LinearAlgebra
using SparseArrays
using AlgebraicMultigrid
using SpecialFunctions: erf

import PyPlot
const plt = PyPlot

include("INC_poisson_3d.jl")

function do_solve( NN::Array{Int64}; stencil_order=3 )
    AA = [-8.0, -8.0, -8.0]
    BB = [ 8.0,  8.0,  8.0]
    
    grid = FD3dGrid( NN, AA, BB )

    # Box dimensions
    Lx = BB[1] - AA[1]
    Ly = BB[2] - AA[2]
    Lz = BB[3] - AA[3]

    # Center of the box
    x0 = AA[1] + Lx/2.0
    y0 = AA[2] + Ly/2.0
    z0 = AA[3] + Lz/2.0

    # Parameters for two gaussian functions
    σ = 0.5
    Npoints = grid.Npoints

    rho = zeros(Float64, Npoints)
    phi = zeros(Float64, Npoints)

    phi_analytic = zeros(Float64, Npoints)
    # Initialization of charge density
    nrmfct = (2*pi*σ^2)^1.5
    #nrmfct = 1.0
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

    println("\nNN = ", NN)
    println("stencil_order = ", stencil_order)
    println("integ rho = ", sum(rho)*dVol)

    lmax = 4
    lmmax = (lmax+1)^2
    Q_lm = zeros(Float64,lmmax)

    calc_multipole_moment!( grid, rho, Q_lm )

    #println("Multipole moment: ")
    #for i in 1:lmmax
    #    @printf("%3d %18.10f\n", i, Q_lm[i])
    #end

    V_x0 = zeros(Float64, grid.Ny, grid.Nz)
    V_xN = zeros(Float64, grid.Ny, grid.Nz)

    V_y0 = zeros(Float64, grid.Nx, grid.Nz)
    V_yN = zeros(Float64, grid.Nx, grid.Nz)

    V_z0 = zeros(Float64, grid.Nx, grid.Ny)
    V_zN = zeros(Float64, grid.Nx, grid.Ny)

    #@time set_bc_isolated!( grid, Q_lm, V_x0, V_xN, V_y0, V_yN, V_z0, V_zN )

    println("Building ∇2")
    ∇2 = build_nabla2_matrix( grid, stencil_order=stencil_order )

    println("Building preconditioner")
    prec = aspreconditioner(ruge_stuben(∇2))

    @printf("Size of ∇2   = %f MiB\n", Base.summarysize(∇2)/(1024*1024))
    @printf("Size of prec = %f MiB\n", Base.summarysize(prec)/(1024*1024))

    print("Solving Poisson equation:\n")
    phi = Poisson_solve_PCG( ∇2, prec, rho, grid, Q_lm, verbose=true )
    #phi = Poisson_solve_PCG( ∇2, prec, rho, 1000 )

    # Calculation of Hartree energy
    Unum = 0.5*sum( rho .* phi ) * dVol
    Uana = 0.5*sum( rho .* phi_analytic ) * dVol
    MAE_pot = sum( abs.( phi .- phi_analytic ) )/Npoints

    rho = reshape(rho, (NN[1],NN[2],NN[3]))
    phi = reshape(phi, (NN[1],NN[2],NN[3]))
    phi_analytic = reshape(phi_analytic, (NN[1],NN[2],NN[3]))

    filepot = open("POT.dat", "w")
    iy = round(Int64,NN[2]/2)
    iz = round(Int64,NN[3]/2)
    println("iy = ", iy)
    println("iz = ", iz)
    for ix in 1:NN[1]
        @printf(filepot, "%18.10f %18.10f %18.10f %18.10f\n",
            grid.x[ix], rho[ix,iy,iz], phi[ix,iy,iz], phi_analytic[ix,iy,iz])
    end
    close(filepot)

    @printf("Numeric  = %18.10f\n", Unum)
    @printf("Uana     = %18.10f\n", Uana)
    @printf("abs diff = %18.10e\n", abs(Unum-Uana))

    @printf("\nMAE pot = %18.10e\n", MAE_pot)

end

using DelimitedFiles
function do_plot(filesave)
    data = readdlm("POT.dat")
    plt.clf()
    plt.plot(data[:,1], data[:,3], label="Numeric")
    plt.plot(data[:,1], data[:,4], label="Analytic")
    plt.grid()
    plt.legend()
    plt.savefig(filesave)
end

function main()
    #
    Ns = [60,60,60]
    #
    do_solve(Ns, stencil_order=3)
    do_plot("IMG_stencil_3.pdf")
    #
    #do_solve(Ns, stencil_order=5)
    #do_plot("IMG_stencil_5.pdf")
    #
    #do_solve(Ns, stencil_order=7)
    #do_plot("IMG_stencil_7.pdf")
    #
    #do_solve(Ns, stencil_order=9)
    #do_plot("IMG_stencil_9.pdf")
end

main()