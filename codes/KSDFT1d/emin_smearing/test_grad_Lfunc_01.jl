push!(LOAD_PATH, "../")

import Random
using Printf
using LinearAlgebra
using Serialization

using KSDFT1d

#include("system_defs_01.jl")
include("system_defs_02.jl")

include("Lfunc.jl")
include("../utilities.jl")
include("gradients_psi_Haux.jl")

# Initialize a Hamiltonian object
Ham = init_Hamiltonian()

hx = Ham.grid.hx
Npoints = Ham.grid.Npoints
Nelectrons = Ham.electrons.Nelectrons
Nstates = Ham.electrons.Nstates
Focc = Ham.electrons.Focc

Random.seed!(1234)
psi = generate_random_wavefunc(Ham)

# Prepare Haux
Haux = psi' * (Ham*psi) * hx # Hsub, subspace Hamiltonian
# Using diagonal Haux
#ebands1 = sort(randn(Nstates))
#Haux = diagm(0 => ebands1)
transform_psi_Haux!(psi, Haux)


# Evaluate total energy by calling Lfunc
E1 = calc_Lfunc_Haux!(Ham, psi, Haux)

g = zeros(Npoints,Nstates)
Hsub = zeros(Nstates,Nstates)
g_Haux = zeros(Nstates,Nstates)
Kg_Haux = zeros(Nstates,Nstates)

# Evaluate gradients
calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)

display(g_Haux)
display(Kg_Haux)