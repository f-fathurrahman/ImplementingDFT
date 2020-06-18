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

function main(N)

    @assert iseven(N)
    
    Random.seed!(1234)
    
    #Ham = create_Ham_H(N, grid_type=:sinc)
    #Ham = create_Ham_LiH(40, grid_type=:sinc)

    #Ham = create_Ham_Ne(N)
    Ham = create_Ham_H2O(N, grid_type=:sinc)

    println(Ham.grid)

    @printf("sizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)

    println("sum V_Ps_loc = ", sum(Ham.V_Ps_loc))

    Nbasis = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi,dVol)

    KS_solve_Emin_PCG!(Ham, psi)

end

main(54)
#for N in [10, 20, 30, 40, 50]
#    main(N)
#end