push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

const DIR_PSP = "../pseudopotentials/pade_gth/"
const DIR_STRUCTURES = "../structures"

include("create_Ham.jl")

include("KS_solve_Emin_PCG.jl")

function main()
    
    Random.seed!(1234)

    #Ham = create_Ham_H2O(51)
    
    #Ham = create_Ham_H(41)
    #Ham = create_Ham_H(40, grid_type=:FD)
    
    #Ham = create_Ham_Ne(41)
    Ham = create_Ham_Ar(41)
    
    #Ham = create_Ham_H2O(40, grid_type=:FD)
    
    #Ham = create_Ham_LiH(40, grid_type=:sinc)

    #Ham = create_Ham_CH4(50, grid_type=:sinc)
    
    #Ham = create_Ham_SiH4(30, grid_type=:sinc)

    #Ham = create_Ham_HCl(30, grid_type=:sinc)

    println(Ham.grid)

    @printf("sizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)

    Nbasis = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi,dVol)

    KS_solve_Emin_PCG!(Ham, psi)

end

@time main()