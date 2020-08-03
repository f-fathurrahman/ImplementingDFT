push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

const DIR_PSP = "../pseudopotentials/pade_gth/"
const DIR_STRUCTURES = "../structures"

include("KS_solve_SCF.jl")
include("KS_solve_Emin_PCG.jl")
include("create_Ham_periodic.jl")

function main()
    
    #Ham = create_Ham_H_periodic(40, grid_type=:FD)
    #Ham = create_Ham_Ne_periodic(60)
    #Ham = create_Ham_LiH_periodic(40)
    #Ham = create_Ham_LiH_periodic_v2(40)
    Ham = create_Ham_CH4_periodic(40, grid_type=:FD)
    #Ham = create_Ham_CH4_periodic(41, grid_type=:LF)

    println(Ham.atoms)
    println(Ham.grid)

    check_betaNL_norm(Ham.grid, Ham.pspotNL)

    @printf("sizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)

    Nbasis = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi,dVol)

    #KS_solve_SCF!(Ham, psi, diag_func=diag_LOBPCG!)
    KS_solve_Emin_PCG!(Ham, psi)
end

@time main()