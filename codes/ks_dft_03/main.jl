push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using KSDFT03Module

const DIR_PSP = "../pseudopotentials/pade_gth/"
const DIR_STRUCTURES = "../structures"

include("smearing.jl")
include("occupations.jl")
include("create_Ham.jl")
include("mix_adaptive.jl")
include("gen_gaussian_density.jl")
include("KS_solve_SCF.jl")
include("KS_solve_SCF_potmix.jl")
include("KS_solve_SCF_NLsolve.jl")

function main( Ham::Hamiltonian; use_smearing=false )
    
    Random.seed!(1234)

    println(Ham.atoms)
    println(Ham.grid)
    check_betaNL_norm(Ham.grid, Ham.pspotNL)

    @printf("\nsizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)

    Nbasis = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi,dVol)

    KS_solve_SCF!(Ham, psi, betamix=0.25, use_smearing=use_smearing,
        guess_density=:random)

    #KS_solve_SCF_potmix!(Ham, psi, betamix=0.25, use_smearing=use_smearing,
    #    guess_density=:random)

    #KS_solve_SCF_NLsolve!(Ham, psi, betamix=0.25, use_smearing=use_smearing,
    #    guess_density=:random)

end

#@time main( create_Ham_Al_atom(40, grid_type=:FD), use_smearing=true )
#@time main( create_Ham_Al_atom(40, grid_type=:LF), use_smearing=true )
#@time main( create_Ham_LiH(40, grid_type=:FD), use_smearing=false )
#@time main( create_Ham_CO(40, grid_type=:FD), use_smearing=false )
#@time main( create_Ham_O2(40, grid_type=:FD), use_smearing=false )
@time main( create_Ham_O2(40, grid_type=:FD, Nstates_extra=2), use_smearing=true )