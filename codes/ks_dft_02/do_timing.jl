push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions
using BenchmarkTools

using KSDFT02Module

const DIR_PSP = "../pseudopotentials/pade_gth/"
const DIR_STRUCTURES = "../structures"

include("create_Ham.jl")

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
    @btime Hpsi = op_H($Ham, $psi)

    println("\nop_K")
    @btime Hpsi = -0.5*$Ham.∇2 * $psi

    println("\nop_V_Ps_nloc")
    @btime Vnlpsi = op_V_Ps_nloc($Ham, $psi)

    #Vnlpsi = zeros(Float64,Nbasis,Nstates)
    #@btime op_V_Ps_nloc!($Ham, $psi, $Vnlpsi)    

    #println("\nbetaNL psi")
    #@btime betaNL_psi = $psi' * $Ham.pspotNL.betaNL *$dVol
    #@time betaNL_psi = psi' * Ham.pspotNL.betaNL *dVol

end
timing_op_H(create_Ham_Al2, 40, :FD)
timing_op_H(create_Ham_Al2, 40, :LF)
