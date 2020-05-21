push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

function main()
    
    atoms = Atoms( xyz_file="NH3.xyz" )
    println(atoms)

    pspfiles = ["H-q1.gth", "N-q5.gth"]

    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [65, 65, 65]

    grid = FD3dGrid( NN, AA, BB )

    println("hx = ", grid.hx)
    println("hy = ", grid.hy)
    println("hz = ", grid.hz)
    println("dVol = ", grid.dVol)
    println(grid.hx*grid.hy*grid.hz)

    Ham = Hamiltonian( atoms, pspfiles, grid )

end

main()