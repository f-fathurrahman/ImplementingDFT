push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

function main()
    atoms = Atoms( xyz_file="NH3.xyz" )
    println(atoms)
end

main()