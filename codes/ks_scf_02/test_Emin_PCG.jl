push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

const DIR_PSP = "../pseudopotentials/pade_gth/"

function create_Ham_H2O()
    atoms = Atoms( xyz_file="H2O.xyz" )
    pspfiles = [ joinpath(DIR_PSP,"O-q6.gth"),
                 joinpath(DIR_PSP, "H-q1.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [41, 41, 41]
    grid = FD3dGrid( NN, AA, BB )
    Ham = Hamiltonian( atoms, pspfiles, grid )
end

function create_Ham_H_atom()
    atoms = Atoms( xyz_string=
        """
        1

        H  0.0  0.0  0.0
        """ )
    pspfiles = ["H-q1.gth"]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [41, 41, 41]
    grid = FD3dGrid( NN, AA, BB )
    return Hamiltonian( atoms, pspfiles, grid )
end

include("KS_solve_Emin_PCG.jl")

function main()
    
    Ham = create_Ham_H2O()
    #Ham = create_Ham_H_atom()

    Nbasis = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi,dVol)

    KS_solve_Emin_PCG!(Ham, psi)

end

main()