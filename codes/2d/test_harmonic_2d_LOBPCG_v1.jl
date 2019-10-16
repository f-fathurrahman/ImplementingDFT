using Printf
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using IncompleteLU

include("init_FD1d_grid.jl")
include("FD2dGrid.jl")
include("build_nabla2_matrix.jl")
include("supporting_functions.jl")

function pot_harmonic( fdgrid::FD2dGrid; ω=1.0 )
    Npoints = fdgrid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = fdgrid.r[1,i]
        y = fdgrid.r[2,i]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 )
    end
    return Vpot
end

function main()
    Nx = 50
    Ny = 50
    fdgrid = FD2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )

    ∇2 = build_nabla2_matrix( fdgrid, func_1d=build_D2_matrix_7pt )

    prec = ilu(-0.5*∇2)

    Vpot = pot_harmonic( fdgrid )
    
    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    # solve for 5 lowest (using `false`) eigenvalues
    res = lobpcg( Ham, false, 5, P=prec )
    #res = lobpcg( Ham, false, 5 )

    println(fieldnames(typeof(res)))

    println(res)
    println(res.converged)
    println("Eigenvalues:")
    for i = 1:5
        @printf("%3d %18.10f %18.10e\n", i, res.λ[i], res.residual_norms[i])
    end
end

main()



