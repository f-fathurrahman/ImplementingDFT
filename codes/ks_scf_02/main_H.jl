push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

const DIR_PSP = "../pseudopotentials/pade_gth/"

include("KS_solve_SCF.jl")

function main()
    
    atoms = Atoms( xyz_string=
        """
        1

        H  0.0  0.0  0.0
        """ )
    println(atoms)

    pspfiles = [joinpath(DIR_PSP,"H-q1.gth")]

    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [41, 41, 41]

    grid = FD3dGrid( NN, AA, BB )
    #grid = LF3dGrid( NN, AA, BB )

    Ham = Hamiltonian( atoms, pspfiles, grid )

    Nbasis = prod(NN)
    Nstates = Ham.electrons.Nstates
    dVol = grid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi)
    psi = psi/sqrt(dVol)

    KS_solve_SCF!( Ham, psi )

end

main()