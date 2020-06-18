using Printf
using LinearAlgebra
using SparseArrays
using Random
using Serialization: serialize

using AlgebraicMultigrid
using SpecialFunctions: erf

include("INC_sch_3d.jl")

include("../common/constants.jl")
include("../common/Atoms.jl")
include("../common/potential_H_atom.jl")

include("../common/ILU0Preconditioner.jl")
include("../common/NoPreconditioner.jl")

function main()

    Random.seed!(1234)

    AA = [-8.0, -8.0, -8.0]
    BB = [ 8.0,  8.0,  8.0]
    NN = [50, 50, 50]

    grid = FD3dGrid( NN, AA, BB )

    atoms = Atoms( xyz_string=
        """
        2

        H   0.75  0.0  0.0
        H  -0.75  0.0  0.0
        """, in_bohr=true)

    V_Ps_loc = pot_Hps_HGH(atoms, grid)

    ∇2 = build_nabla2_matrix( grid )

    Ham = -0.5*∇2 + spdiagm( 0 => V_Ps_loc )

    println("Building preconditioner")
    #prec = aspreconditioner(ruge_stuben(-0.5*∇2))
    #prec = aspreconditioner(ruge_stuben(Ham))
    prec = ILU0Preconditioner(Ham)
    #prec = ILU0Preconditioner(-0.5*∇2)
    println("Done building preconditioner")

    Nstates = 5  # only choose the lowest lying state
    Npoints = grid.Npoints
    psi = rand(Float64, Npoints, Nstates)
    dVol = grid.dVol
    ortho_sqrt!(psi,dVol)

    evals = diag_LOBPCG!( Ham, psi, prec, verbose=true )
    #evals = diag_Emin_PCG!( Ham, X, prec, verbose=true )

    psi = psi/sqrt(dVol) # renormalize
    serialize("wavefunc.data", psi)
end

main()
