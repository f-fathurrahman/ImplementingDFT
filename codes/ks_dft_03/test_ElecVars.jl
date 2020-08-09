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
include("gen_gaussian_density.jl")

include("ElecVars.jl")
include("emin_smearing.jl")

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
    kT = 0.01

    evars = ElecVars(Ham, psi)
    println(evars)

    Ham.electrons.eorbs[:] = evars.Hsub_eigs[:]

    g = ElecGradient(Ham)
    Kg = ElecGradient(Ham)

    calc_energies_grad!( Ham, evars, g, Kg, kT )

    println("evars.Hsub")
    display(evars.Hsub); println()
    println("g.Haux")
    display(g.Haux); println()
    println("Kg.Haux")
    display(Kg.Haux); println()

end

main( create_Ham_Al_atom(40, grid_type=:FD), use_smearing=true )
#main( create_Ham_C_atom(40, grid_type=:FD), use_smearing=true )