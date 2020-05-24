push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

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

function main(;N=41)
    
    #Ham = create_Ham_NH3(N)
    Ham = create_Ham_CH4(N)

    dVol = Ham.grid.dVol
    NbetaNL = Ham.pspotNL.NbetaNL
    betaNL = Ham.pspotNL.betaNL
    for ibeta in 1:NbetaNL
        @printf("%3d %18.10f\n", ibeta, sum(betaNL[:,ibeta].*betaNL[:,ibeta])*dVol)
    end

end

main(N=31)
main(N=41)
main(N=51)