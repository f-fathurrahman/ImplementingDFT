using Printf
using LinearAlgebra
using SparseArrays
using SpecialFunctions

using AlgebraicMultigrid
using Random
using Serialization

include("INC_sch_3d.jl")

# Caution: make sure to avoid the singularity
function pot_H_atom( grid::FD3dGrid; r0=(0.0, 0.0, 0.0) )
    Npoints = grid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        dx = grid.r[1,i] - r0[1]
        dy = grid.r[2,i] - r0[2]
        dz = grid.r[3,i] - r0[3]
        Vpot[i] = -1.0/sqrt(dx^2 + dy^2 + dz^2)
    end
    return Vpot
end

function pot_Hps_GTH( grid::FD3dGrid; r0=(0.0, 0.0, 0.0) )
    Npoints = grid.Npoints
    Vpot = zeros( Float64, Npoints )

    # Parameters
    Zval = 1
    rloc = 0.2
    C1 = -4.0663326
    C2 = 0.6678322

    # TODO Add journal reference
    for ip = 1:Npoints
        dx2 = ( grid.r[1,ip] - r0[1] )^2
        dy2 = ( grid.r[2,ip] - r0[2] )^2
        dz2 = ( grid.r[3,ip] - r0[3] )^2
        r = sqrt(dx2 + dy2 + dz2)
        if r < eps()
            Vpot[ip] = -2*Zval/(sqrt(2*pi)*rloc) + C1
        else
            rrloc = r/rloc
            Vpot[ip] = -Zval/r * erf( r/(sqrt(2.0)*rloc) ) +
                     (C1 + C2*rrloc^2)*exp(-0.5*(rrloc)^2)
        end
    end
    return Vpot
end

function main()

    Random.seed!(1234)

    Nx = 50
    Ny = 50
    Nz = 50
    A = 6.0
    grid = FD3dGrid( (-A,A), Nx, (-A,A), Ny, (-A,A), Nz )

    ∇2 = build_nabla2_matrix( grid, stencil_order=11 )

    Vpot = pot_H_atom( grid )
    #Vpot = pot_Hps_GTH( grid )

    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    println("Building preconditioner")
    prec = aspreconditioner(ruge_stuben(Ham))
    println("Done building preconditioner")

    Nstates = 1  # only choose the lowest lying state
    Npoints = Nx*Ny*Nz
    X = rand(Float64, Npoints, Nstates)
    ortho_sqrt!(X)

    evals = diag_LOBPCG!( Ham, X, prec, verbose=true )

    @printf("\n\nEigenvalues\n")
    for i in 1:Nstates
        @printf("%5d %18.10f\n", i, evals[i])
    end
end

main()



