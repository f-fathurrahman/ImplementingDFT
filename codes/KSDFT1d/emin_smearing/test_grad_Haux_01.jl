push!(LOAD_PATH, "../")

import Random
using Printf
using LinearAlgebra
using Serialization

using KSDFT1d

include("system_defs_01.jl")
include("Lfunc.jl")
include("../utilities.jl")
include("gradients_psi_Haux.jl")

Ham = init_Hamiltonian()

hx = Ham.grid.hx
Npoints = Ham.grid.Npoints
Nelectrons = Ham.electrons.Nelectrons
Nstates = Ham.electrons.Nstates
Focc = Ham.electrons.Focc

Random.seed!(1234)
psi = generate_random_wavefunc(Ham)

Haux = psi' * (Ham*psi) * hx # Hsub, subspace Hamiltonian

E1 = calc_Lfunc_Haux!(Ham, psi, Haux)

g = zeros(Npoints,Nstates)
Hsub = zeros(Nstates,Nstates)
calc_grad!(Ham, psi, g, Hsub)

g_Haux = zeros(Nstates,Nstates)
Kg_Haux = zeros(Nstates,Nstates)
calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)