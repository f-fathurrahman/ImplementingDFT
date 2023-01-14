push!(LOAD_PATH, "../")

import Random
using Printf
using LinearAlgebra
using Serialization

using KSDFT1d

include("system_defs_01.jl")
include("Lfunc.jl")
include("../utilities.jl")

Ham = init_Hamiltonian()

hx = Ham.grid.hx
Npoints = Ham.grid.Npoints
Nelectrons = Ham.electrons.Nelectrons
Nstates = Ham.electrons.Nstates
Focc = Ham.electrons.Focc

Random.seed!(1234)
psi = generate_random_wavefunc(Ham)

# Generate Haux
Haux = Hermitian(rand(Nstates,Nstates))

_, Urot = eigen(Haux)
Etot1 = calc_Lfunc_Haux!(Ham, psi, Haux)

println("\nCall again: (by simultaneously transforming psi and Haux)")
Etot2 = calc_Lfunc_Haux!(Ham, psi*Urot, Urot' * Haux * Urot)