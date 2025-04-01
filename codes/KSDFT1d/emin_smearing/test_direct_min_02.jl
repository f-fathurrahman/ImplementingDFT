import Random
using Printf
using LinearAlgebra
using Serialization

using KSDFT1d

include("system_defs_01.jl")
include("KohnSham_Etotal.jl")

function prec_invK(Ham::Hamiltonian1d, v)
    return inv(Ham.Kmat)*v
end

include("../utilities.jl")

#=
function main()

    Ham = init_Hamiltonian()

    hx = Ham.grid.hx
    Npoints = Ham.grid.Npoints
    Nelectrons = Ham.electrons.Nelectrons
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc

    Kmat = Ham.Kmat
    Vtot = Ham.potentials.Total
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    rhoe = Ham.rhoe

    Nspin = 1
    Hmat = zeros(Float64, Npoints, Npoints)
    epsxc = zeros(Float64, Npoints)
    rhoe_new = zeros(Float64, Npoints, Nspin)

    Etot = Inf
    Etot_old = Etot
    E_NN = calc_E_NN(Ham.atoms)

    Ekin = 0.0
    Ehartree = 0.0
    Eion = 0.0
    Exc = 0.0

    betamix = 0.5
    mixer = BroydenMixer(rhoe, betamix, mixdim=8)

    Focc = Ham.electrons.Focc
    use_smearing = true
    kT = 0.1*eV2Ha #0.1 # 0.3*eV2Ha
    mTS = 0.0

    println("kT = ", kT)

    evals = Ham.electrons.ebands


end

main()
=#