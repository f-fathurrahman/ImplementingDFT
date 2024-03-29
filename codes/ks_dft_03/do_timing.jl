push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions
using BenchmarkTools

using KSDFT03Module

const DIR_PSP = "../pseudopotentials/pade_gth/"
const DIR_STRUCTURES = "../structures"

include("../ks_dft_02/create_Ham.jl") # for testing only

function timing_create_Ham(Ham_func, N, grid_type)
    @time Ham = Ham_func(N, grid_type=grid_type)
    @time Ham = Ham_func(N, grid_type=grid_type)
    @printf("\nsizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)
end
#timing_create_Ham(create_Ham_Al2, 40, :FD)


function timing_op_H(Ham_func, N, grid_type)
    Ham = Ham_func(N, grid_type=grid_type)
    Nbasis = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol
    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi,dVol)

    println("\nop_H")
    @time Hpsi = op_H(Ham, psi)
    @time Hpsi = op_H(Ham, psi)

    println("\nop_K")
    @time Hpsi = -0.5*Ham.Laplacian * psi
    @time Hpsi = -0.5*Ham.Laplacian * psi

    println("\nop_V_Ps_nloc")
    @btime Vnlpsi = op_V_Ps_nloc($Ham, $psi)

    #println("\nbetaNL psi")
    #@time betaNL_psi = psi' * Ham.pspotNL.betaNL *dVol
    #@time betaNL_psi = psi' * Ham.pspotNL.betaNL *dVol

end
timing_op_H(create_Ham_Al2, 40, :FD)
