push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

const DIR_PSP = "../pseudopotentials/pade_gth/"

include("create_Ham.jl")

include("KS_solve_Emin_PCG.jl")

function main()
    
    #Ham = create_Ham_H2O(41)
    Ham = create_Ham_H(41)
    #Ham = create_Ham_Ne(41)

    Nbasis = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi,dVol)

    KS_solve_Emin_PCG!(Ham, psi)

end

main()