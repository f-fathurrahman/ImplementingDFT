push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using KSDFT02Module

const DIR_PSP = "../pseudopotentials/pade_gth/"
const DIR_STRUCTURES = "../structures"

include("create_Ham.jl")
include("KS_solve_SCF.jl")
include("KS_solve_SCF_potmix.jl")
include("KS_solve_Emin_PCG.jl")
include("KS_solve_TRDCM.jl")

function main(Ham_func, N, grid_type)
    
    Random.seed!(1234)

    Ham = Ham_func(N, grid_type=grid_type)

    #for isp in 1:Ham.atoms.Nspecies
    #    println(Ham.pspots[isp])
    #end

    @printf("\nsizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)

    Nbasis = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol
    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi,dVol)

    println("Some psi")
    println(psi[1,1])
    println(psi[2,1])
    println(psi[3,1])

    Vpsnl = op_V_Ps_nloc(Ham, psi)
    println("dot psi,Vpsnl = ", dot(psi,Vpsnl))

    @time Epsnl = calc_E_Ps_nloc(Ham, psi)
    @time Epsnl = calc_E_Ps_nloc(Ham, psi)
    println("Epsnl = ", Epsnl)
end

#@time main(create_Ham_CO, 40, :FD)
#@time main(create_Ham_H2O, 40, :FD)
@time main(create_Ham_LiH, 40, :FD)
#@time main(create_Ham_NH3, 40, :FD)
#@time main(create_Ham_CH4, 40, :FD)
#@time main(create_Ham_Al2, 40, :FD)
#@time main(create_Ham_Ni2, 40, :FD)
