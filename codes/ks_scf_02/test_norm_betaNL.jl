push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

const DIR_PSP = "../pseudopotentials/pade_gth/"

function create_Ham_NH3( N::Int64 )
    atoms = Atoms( xyz_file="NH3.xyz" )
    pspfiles = ["N-q5.gth", "H-q1.gth"]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N, N, N]
    grid = FD3dGrid( NN, AA, BB )

    return Hamiltonian( atoms, pspfiles, grid )
end

function create_Ham_CH4( N::Int64 )
    atoms = Atoms( xyz_file="CH4.xyz" )
    pspfiles = ["C-q4.gth", "H-q1.gth"]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    grid = FD3dGrid( NN, AA, BB )
    
    return Hamiltonian( atoms, pspfiles, grid )
end

function create_Ham_Ar( N::Int64 )
    atoms = Atoms( xyz_string=
        """
        1

        Ar  0.0  0.0  0.0
        """ )
    pspfiles = [ joinpath(DIR_PSP, "Ar-q8.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    grid = FD3dGrid( NN, AA, BB )

    return Hamiltonian( atoms, pspfiles, grid )
end

function create_Ham_Ne( N::Int64 )
    atoms = Atoms( xyz_string=
        """
        1

        Ne  0.0  0.0  0.0
        """ )
    pspfiles = ["Ne-q8.gth"]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    grid = FD3dGrid( NN, AA, BB )

    return Hamiltonian( atoms, pspfiles, grid )
end

function create_Ham_SiH4(N::Int64)
    atoms = Atoms( xyz_file="SiH4.xyz" )
    pspfiles = ["Si-q4.gth", "H-q1.gth"]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    grid = FD3dGrid( NN, AA, BB )

    return Hamiltonian( atoms, pspfiles, grid )
end

function main(;N=41)
    
    #Ham = create_Ham_NH3(N)
    #Ham = create_Ham_CH4(N)
    Ham = create_Ham_Ne(N)
    #Ham = create_Ham_Ar(N)
    #Ham = create_Ham_SiH4(N)

    dVol = Ham.grid.dVol
    NbetaNL = Ham.pspotNL.NbetaNL
    betaNL = Ham.pspotNL.betaNL
    println("Test normalization")
    for ibeta in 1:NbetaNL
        @printf("%3d %18.10f\n", ibeta, sum(betaNL[:,ibeta].*betaNL[:,ibeta])*dVol)
    end

end

main(N=31)
main(N=41)
main(N=51)
main(N=61)