push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

const DIR_PSP = "../pseudopotentials/pade_gth/"
const DIR_STRUCTURES = "../structures"

include("smearing.jl")
include("occupations.jl")
include("create_Ham.jl")
include("mix_adaptive.jl")
include("mix_linear.jl")
include("mix_rpulay.jl")
include("gen_gaussian_density.jl")
include("KS_solve_SCF.jl")

function main( Ham::Hamiltonian; use_smearing=false )
    
    Random.seed!(1234)

    println(Ham.atoms)
    println(Ham.grid)
    check_betaNL_norm(Ham.grid, Ham.pspotNL)

    @printf("\nsizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)

    Nbasis = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol

    Nspin = Ham.Nspin
    psis = Vector{Matrix{Float64}}(undef,Nspin)
    for i in 1:Nspin
        psis[i] = rand(Float64, Nbasis, Nstates)
        ortho_sqrt!(psis[i], dVol)
    end

    KS_solve_SCF!(Ham, psis, betamix=0.5, use_smearing=use_smearing,
        guess_density=:random)
end

#@time main( create_Ham_O2(40, grid_type=:FD, Nstates_extra=2), use_smearing=true )

#@time main(
#    create_Ham_Al_atom(40, grid_type=:FD, Nstates_extra=2),
#    use_smearing=false
#)

#@time main(
#    create_Ham_Al_atom(40, grid_type=:FD, Nstates_extra=0),
#    use_smearing=false
#)

@time main(
    create_Ham_Al2(40, grid_type=:FD, Nstates_extra=0),
    use_smearing=false
)

#@time main(
#    create_Ham_C_atom(40, grid_type=:FD, Nstates_extra=0),
#    use_smearing=false
#)

#@time main(
#    create_Ham_Al_atom(40, grid_type=:FD, Nstates_extra=4),
#    use_smearing=true, manual_Focc=true
#)