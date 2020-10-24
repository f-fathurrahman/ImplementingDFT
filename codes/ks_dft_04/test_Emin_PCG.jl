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
include("KS_solve_SCF_potmix.jl")
#include("KS_solve_Emin_PCG.jl")
include("KS_solve_Emin_PCG_v2.jl")

function main( N::Int64; grid_type=:FD, use_smearing=false, kT=1.e-3 )

    Random.seed!(1234)

    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES,"O2.xyz") )
    pspfiles = [ joinpath(DIR_PSP,"O-q6.gth") ]

    #atoms = Atoms( xyz_string=
    #    """
    #    2
    #    
    #    C   0.575  0.0  0.0
    #    O  -0.575  0.0  0.0
    #    """) # coordinates are in angstrom
    #pspfiles = [ joinpath(DIR_PSP,"C-q4.gth"),
    #             joinpath(DIR_PSP,"O-q6.gth") ]

    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    if (grid_type == :sinc) || (grid_type == :LF)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB )
    end
    #Ham = Hamiltonian( atoms, pspfiles, grid, N_unpaired=0, Nspin=1 )
    Ham = Hamiltonian( atoms, pspfiles, grid, N_unpaired=2, Nspin=1 )
    #Ham = Hamiltonian( atoms, pspfiles, grid, N_unpaired=0, Nspin=2 )

    #Ham = Hamiltonian( atoms, pspfiles, grid, N_unpaired=0, Nspin=2, Nstates_extra=2 )

    #Ham = Hamiltonian( atoms, pspfiles, grid, N_unpaired=2, Nspin=2 )
    #Ham = Hamiltonian( atoms, pspfiles, grid, N_unpaired=0, Nstates_extra=1, Nspin=1 )
    #Ham = Hamiltonian( atoms, pspfiles, grid, Nstates_extra=2, Nspin=1 )

    println(Ham.atoms)
    println(Ham.grid)
    check_betaNL_norm(Ham.grid, Ham.pspotNL)

    @printf("\nsizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)

    Npoints = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol

    Nspin = Ham.Nspin
    psis = Vector{Matrix{Float64}}(undef,Nspin)
    for i in 1:Nspin
        psis[i] = rand(Float64, Npoints, Nstates)
        ortho_sqrt!(psis[i], dVol)
    end

    KS_solve_Emin_PCG!(Ham, psis)
    #KS_solve_SCF!(Ham, psis, use_smearing=true, kT=1e-3, betamix=0.2)
    #KS_solve_SCF!(Ham, psis, betamix=0.1)
    #KS_solve_SCF_potmix!(Ham, psis, use_smearing=true, kT=1e-3, betamix=0.2)

end

@time main(40, grid_type=:FD)
