using Printf
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using IncompleteLU
using AlgebraicMultigrid
using Random
using Serialization

include("INC_sch_3d_LF.jl")

function pot_harmonic( lfgrid::LF3dGrid; ω=1.0 )
    Npoints = lfgrid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = lfgrid.r[1,i]
        y = lfgrid.r[2,i]
        z = lfgrid.r[3,i]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 + z^2 )
    end
    return Vpot
end

function main()

    Random.seed!(1234)

    Nx = 25
    Ny = 25
    Nz = 25
    lfgrid = LF3dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny, (-5.0,5.0), Nz,
        type_x=:sinc, type_y=:sinc, type_z=:sinc )

    ∇2 = build_nabla2_matrix( lfgrid )

    Vpot = pot_harmonic( lfgrid )
    
    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    # may choose between these two
    println("Building preconditioner")
    #prec = ilu(-0.5*∇2)
    #prec = ilu(Ham) # this should result in faster convergence
    prec = aspreconditioner(ruge_stuben(Ham))
    #prec = aspreconditioner(smoothed_aggregation(Ham))
    #prec = aspreconditioner(ruge_stuben(-0.5*∇2))
    #prec = aspreconditioner(smoothed_aggregation(-0.5*∇2))
    println("Done building preconditioner")

    Nstates = 10
    Npoints = Nx*Ny*Nz
    X = rand(Float64, Npoints, Nstates)
    ortho_sqrt!(X)
    
    #evals = diag_Emin_PCG!( Ham, X, prec, verbose=true )
    #evals = diag_davidson!( Ham, X, prec, verbose=true )
    evals = diag_LOBPCG!( Ham, X, prec, verbose=true )

    @printf("\n\nEigenvalues\n")
    for i in 1:Nstates
        @printf("%5d %18.10f\n", i, evals[i])
    end

    #serialize("wavefunc.data", X)

end

main()



