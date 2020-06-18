push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

const DIR_PSP = "../pseudopotentials/pade_gth/"

include("create_Ham.jl")

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