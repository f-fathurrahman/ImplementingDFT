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
include("KS_solve_Emin_PCG.jl")

function main( N::Int64; grid_type=:FD, use_smearing=false, kT=1.e-3 )

    Random.seed!(1234)

    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES,"O2.xyz") )
    pspfiles = [ joinpath(DIR_PSP,"O-q6.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    if (grid_type == :sinc) || (grid_type == :LF)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB )
    end
    Ham = Hamiltonian( atoms, pspfiles, grid, N_unpaired=0, Nspin=1 )

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

    println(Ham.electrons)

    g = Vector{Matrix{Float64}}(undef,Nspin)
    Kg = Vector{Matrix{Float64}}(undef,Nspin)
    Hsub = Vector{Matrix{Float64}}(undef,Nspin)
    for i in 1:Nspin
        g[i] = zeros(Float64, Npoints, Nstates)
        Kg[i] = zeros(Float64, Npoints, Nstates)
        Hsub[i] = zeros(Float64, Nstates, Nstates)
    end

    Rhoe = zeros(Float64, Npoints, Nspin)

    calc_rhoe!( Ham, psis, Rhoe )
    update!( Ham, Rhoe )
    calc_energies!( Ham, psis )

    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin

    for i in 1:Nspin
        Ham.ispin = i
        calc_grad!( Ham, psis[i], g[i], Hsub[i] )
        # apply preconditioner
        Kg[i][:,:] = g[i][:,:]
        for ist in 1:Nstates
            @views ldiv!(Ham.precKin, Kg[i][:,ist])
        end
    end

    println(Ham.energies)

end

@time main(40, grid_type=:FD)
