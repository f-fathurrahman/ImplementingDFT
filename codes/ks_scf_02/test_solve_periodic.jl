push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

const DIR_PSP = "../pseudopotentials/pade_gth/"
const DIR_STRUCTURES = "../structures"

include("KS_solve_SCF.jl")
include("KS_solve_Emin_PCG.jl")

function create_Ham_H( N::Int64; grid_type=:FD, pbc=(false,false,false) )
    atoms = Atoms( xyz_string=
        """
        1

        H  8.0  8.0  8.0
        """, in_bohr=true, pbc=pbc, LatVecs=16.0*diagm(ones(3)) )
    pspfiles = [joinpath(DIR_PSP,"H-q1.gth")]
    AA = zeros(3)
    BB = 16.0*ones(3)
    NN = [N,N,N]
    if grid_type == :sinc
        @assert pbc == (false,false,false)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB, pbc=pbc )
    end
    return Hamiltonian( atoms, pspfiles, grid )
end


function main()
    
    Ham = create_Ham_H(64, pbc=(true,true,true))

    println(Ham.atoms)

    @printf("sizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)

    println(Ham.grid)

    Nbasis = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi,dVol)

    #KS_solve_SCF!(Ham, psi, diag_func=diag_LOBPCG!)
    KS_solve_Emin_PCG!(Ham, psi)
end

main()