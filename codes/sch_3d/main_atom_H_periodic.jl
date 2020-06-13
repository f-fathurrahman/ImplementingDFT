using Printf
using LinearAlgebra
using SparseArrays
using FFTW

using AlgebraicMultigrid
using Random

include("INC_sch_3d.jl")

include("../common/Atoms.jl")
include("../common/PsPot_GTH.jl")
include("../common/GVectors.jl")
include("../common/calc_strfact.jl")
include("../common/XSF_utils.jl")

function init_V_Ps_loc_G( atpos, grid, gvec )

    Npoints = grid.Npoints
    CellVolume = grid.Lx * grid.Ly * grid.Lz  # FIXME: orthogonal LatVecs
    Nx = grid.Nx
    Ny = grid.Ny
    Nz = grid.Nz
    G2 = gvec.G2
    Ng = length(G2)

    V_Ps_loc = zeros(Float64, Npoints)
    Vg = zeros(ComplexF64, Npoints)
    
    strf = calc_strfact( atpos, 1, [1], gvec.G )

    psp = PsPot_GTH("../pseudopotentials/pade_gth/H-q1.gth")
    isp = 1
    for ig = 1:Ng
        Vg[ig] = strf[ig,isp] * eval_Vloc_G( psp, G2[ig] )
    end
    #
    ctmp = reshape(Vg, (Nx,Ny,Nz) )
    ifft!(ctmp)
    ctmp = reshape(ctmp, Npoints)
    V_Ps_loc[:] = V_Ps_loc[:] + real( ctmp ) * Npoints / CellVolume

    return V_Ps_loc
end

function main(N)

    Random.seed!(1234)

    Nx = N
    Ny = N
    Nz = N
    L = 16.0
    grid = FD3dGrid( (0.0,L), Nx, (0.0,L), Ny, (0.0,L), Nz, pbc=(true,true,true) )

    #println(grid)

    ∇2 = build_nabla2_matrix( grid )

    atpos = zeros(3,1)
    atpos[:,1] = [0.0, 0.0, 0.0]

    gvec = GVectors(grid)
    Vpot = init_V_Ps_loc_G( atpos, grid, gvec )

    #filnam = "OUT_Vpot_atom_H.xsf"
    #LatVecs = [L 0.0 0.0; 0.0 L 0.0; 0.0 0.0 L]
    #Ns = (Nx, Ny, Nz)
    #write_xsf(filnam, LatVecs, atpos)
    #write_xsf_data3d_crystal(filnam, Ns, LatVecs, Vpot)

    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    #println("Building preconditioner")
    prec = aspreconditioner(ruge_stuben(Ham))
    #println("Done building preconditioner")

    Nstates = 1  # only choose the lowest lying state
    Npoints = Nx*Ny*Nz
    X = rand(Float64, Npoints, Nstates)
    ortho_sqrt!(X)

    evals = diag_LOBPCG!( Ham, X, prec, verbose=false )

    @printf("\n\nEigenvalues (N=%d) h = %18.10f\n", N, grid.hx)
    for i in 1:Nstates
        @printf("%5d %18.10f\n", i, evals[i])
    end
end

for N in [31,35,41,45,51]
    main(N)
end




