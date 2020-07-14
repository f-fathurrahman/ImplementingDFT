using Printf
using LinearAlgebra
using SparseArrays
using Random

include("INC_qdot_2d.jl")

function pot_harmonic( grid::FD2dGrid; ω=1.0 )
    Npoints = grid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = grid.r[1,i]
        y = grid.r[2,i]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 )
    end
    return Vpot
end

function do_run( N::Int64; Nstates=1 )

    Random.seed!(1234)

    Nx = N
    Ny = N
    L = 50.0
    grid = FD2dGrid( (-L/2,L/2), Nx, (-L/2,L/2), Ny )

    ∇2 = build_nabla2_matrix( grid )

    Vpot = pot_harmonic( grid, ω=0.22 )
    
    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    prec = ILU0Preconditioner(Ham)

    dVol = grid.dVol
    Npoints = Nx*Ny
    X = rand(Float64, Npoints, Nstates)
    ortho_sqrt!(X, dVol)

    #evals = diag_Emin_PCG!( Ham, X, prec )
    evals = diag_LOBPCG!( Ham, X, prec )

    X = X/sqrt(grid.dVol) # renormalize

    @printf("%5d ", N)
    for i in 1:Nstates
        @printf(" %18.10f", evals[i])
    end
    @printf("\n")

end

for N in range(10, stop=60, step=5)
    do_run(N, Nstates=3)
end


