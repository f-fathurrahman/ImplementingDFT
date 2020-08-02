push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

const DIR_PSP = "../pseudopotentials/pade_gth/"
const DIR_STRUCTURES = "../structures"

include("create_Ham.jl")

function main( create_Ham_func; N=41)    
    Ham = create_Ham_func(N)
    check_betaNL_norm(Ham.grid, Ham.pspotNL)
end

for N in range(21,stop=61,step=10)
    main(create_Ham_CH4, N=N)
end
