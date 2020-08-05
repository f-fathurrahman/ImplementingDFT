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

function main(Ham_func, N, grid_type)
    
    Random.seed!(1234)

    Ham = Ham_func(N, grid_type=grid_type)

    println(Ham.atoms)
    println(Ham.grid)
    check_betaNL_norm(Ham.grid, Ham.pspotNL)

    @printf("\nsizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)

    Nbasis = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi,dVol)

    KS_solve_Emin_PCG!(Ham, psi)

end

#for N in [21, 25, 31, 35, 41, 45, 51, 55, 61]
#    main(create_Ham_H, N, :FD)
#end

for N in [20, 25, 30, 35, 40, 45, 50, 55, 60]
    main(create_Ham_CO, N, :LF)
end

#for N in [20, 26, 30, 36, 40, 46, 50, 56, 60]
#    main(create_Ham_H2, N, :FD)
#end
