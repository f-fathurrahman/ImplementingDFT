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
include("gen_gaussian_density.jl")
include("mix_adaptive.jl")
include("KS_solve_SCF.jl")
include("KS_solve_SCF_NLsolve.jl")

include("ElecVars.jl")
include("emin_smearing.jl")
include("linmin_grad.jl")
include("KS_solve_Emin_SD_Haux.jl")

function main( Ham::Hamiltonian; use_smearing=false )
    
    Random.seed!(1234)

    Nbasis = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi,dVol)
    kT = 0.01

    evars = ElecVars(Ham, psi)
    println(evars)
    #Ham.electrons.eorbs[:] = [-6.0, -5.0, -4.0, -3.0, -2.0, -1.0] #evars.Hsub_eigs[:]

    println("Initial eorbs is set to Hsub_eigs")
    Ham.electrons.eorbs[:] = evars.Hsub_eigs[:]

    g = ElecGradient(Ham)
    Kg = ElecGradient(Ham)
    subrot = SubspaceRotations(Nstates)
    
    Ham.energies.NN = calc_E_NN( Ham.atoms, Ham.pspots )
    
    Etot = compute!( Ham, evars, g, Kg, kT, subrot )
    ss = dot(evars.psi, g.psi)*dVol
    @printf("dot evars.psi and g.psi     = %18.10f\n", ss)

    ss = dot(diagm(0 => Ham.electrons.eorbs), g.Haux)
    @printf("dot diagm(eorbs) and g.Haux = %18.10f\n", ss)

    # Doing step
    d = deepcopy(g)
    d.psi[:] = -g.psi[:]
    d.Haux[:]  = -g.Haux[:]
    α = 0.0
    α_Haux = 0.1

    do_step!( Ham, α, α_Haux, evars, d, subrot )

    Etot_new = compute!( Ham, evars, g, Kg, kT, subrot )    

    ss = dot(evars.psi, g.psi)*dVol
    @printf("dot evars.psi and g.psi     = %18.10f\n", ss)

    ss = dot(diagm(0 => Ham.electrons.eorbs), g.Haux)
    @printf("dot diagm(eorbs) and g.Haux = %18.10f\n", ss)

end

main( create_Ham_Al_atom(40, grid_type=:FD), use_smearing=true )
#main( create_Ham_C_atom(40, grid_type=:FD), use_smearing=true )