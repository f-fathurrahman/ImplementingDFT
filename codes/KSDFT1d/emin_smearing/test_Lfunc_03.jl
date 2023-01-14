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

# Using diagonal Haux
ebands1 = sort(randn(Nstates))
Haux = diagm(0 => ebands1)

display(Haux); println()

Etot1 = calc_Lfunc_Haux!(Ham, psi, Haux)

println("\nCall again:")
Etot2 = calc_Lfunc_ebands!(Ham, psi, reshape(ebands1,(Nstates,1)))