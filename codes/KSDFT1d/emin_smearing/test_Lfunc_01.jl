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
Hsub = psi' * (Ham*psi) * hx

ebands = zeros(Nstates,1) # explicitly declare ebands as 1-column matrix
ebands[:,1], Urot = eigen(Hsub)

display(ebands); println()

Etot1 = calc_Lfunc_ebands!(Ham, psi, ebands)

println("\nUsing rotated psi")
psi2 = psi*Urot
Etot2 = calc_Lfunc_ebands!(Ham, psi*Urot, ebands)

#display(psi2' * psi2 * hx); println()