push!(LOAD_PATH, pwd())

using Printf
using BenchmarkTools
using MyModule

using Random
Random.seed!(1234)

function pot_harmonic( grid; ω=1.0 )
    Npoints = grid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = grid.r[1,i]
        y = grid.r[2,i]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 )
    end
    return Vpot
end

function time_init_grid(N, L::Float64)
    @printf("Timing FD2dGrid, N=%d: ", N)
    @btime grid = FD2dGrid( (-$L/2,$L/2), $N, (-$L/2,$L/2), $N )
end

function time_build_nabla2(N, L::Float64)
    @printf("Timing build_nabla2_matrix, N=%d: ", N)
    grid = FD2dGrid( (-L/2,L/2), N, (-L/2,L/2), N )
    @btime ∇2 = build_nabla2_matrix( $grid )
end

function time_Hamiltonian(N, L::Float64; Nstates=1)
    @printf("Timing build_Hamiltonian, N=%d, Nstates=%d: ", N, Nstates)
    grid = FD2dGrid( (-L/2,L/2), N, (-L/2,L/2), N )
    V_ext = pot_harmonic( grid, ω=0.22 )
    Nelectrons = 2*Nstates
    @btime Ham = Hamiltonian( $grid, $V_ext, Nelectrons=$Nelectrons )
end

function time_op_H(N, L::Float64; Nstates=1)
    grid = FD2dGrid( (-L/2,L/2), N, (-L/2,L/2), N )
    V_ext = pot_harmonic( grid, ω=0.22 )
    #
    Nelectrons = 2*Nstates
    Ham = Hamiltonian( grid, V_ext, Nelectrons=Nelectrons )
    Nbasis = Ham.grid.Npoints
    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi)
    psi = psi/sqrt(grid.dVol)
    #
    @printf("Timing op_H, N=%d, Nstates=%d: ", N, Nstates)
    @btime Hpsi = op_H($Ham, $psi)
    #
    Hpsi = zeros(Float64,Nbasis,Nstates)
    @printf("Timing op_H, N=%d, Nstates=%d (in-place): ", N, Nstates)
    @btime op_H!($Ham, $psi, $Hpsi)
end


for N in [20, 40, 80]
    println()
    time_init_grid(N, 20.0)
    time_build_nabla2(N, 20.0)
    time_Hamiltonian(N, 20.0, Nstates=2)
    time_op_H(N, 20.0, Nstates=2)
end