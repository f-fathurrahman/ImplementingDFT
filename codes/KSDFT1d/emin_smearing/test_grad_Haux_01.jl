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

dx = Ham.grid.dx
Npoints = Ham.grid.Npoints
Nelectrons = Ham.electrons.Nelectrons
Nstates = Ham.electrons.Nstates
Focc = Ham.electrons.Focc

Random.seed!(1234)
psi = generate_random_wavefunc(Ham)

Haux = psi' * (Ham*psi) * dx # Hsub, subspace Hamiltonian

# Evaluate total energy by calling Lfunc
E1 = calc_Lfunc_Haux!(Ham, psi, Haux)

g = zeros(Npoints,Nstates)
Hsub = zeros(Nstates,Nstates)
# Evaluate the gradient for psi
calc_grad!(Ham, psi, g, Hsub)
# The outputs are g and Hsub
# Hsub is byproduct of calculation of g
# Hsub will be used in calculating g_Haux and Kg_Haux

# Evaluate the gradient for Haux
g_Haux = zeros(Nstates,Nstates)
Kg_Haux = zeros(Nstates,Nstates)
calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
# The additional input is Hsub
# The outputs are in g_Haux and Kg_Haux

# FIXME: these should be combined into one function call
