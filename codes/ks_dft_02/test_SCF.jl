push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

const DIR_PSP = "../pseudopotentials/pade_gth/"
const DIR_STRUCTURES = "../structures"

include("create_Ham.jl")

include("KS_solve_SCF.jl")

function main()
    
    #Ham = create_Ham_H2O(50)
    #Ham = create_Ham_H(41)
    #Ham = create_Ham_H2(50)
    #Ham = create_Ham_Ne(41)
    Ham = create_Ham_CH4(41)
    #Ham = create_Ham_LiH(41)

    @printf("sizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)

    println(Ham.grid)

    Nbasis = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi,dVol)

    KS_solve_SCF!(Ham, psi, diag_func=diag_LOBPCG!)

end

@time main()