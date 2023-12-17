push!(LOAD_PATH, "../")

import Random
using Printf
using LinearAlgebra
using Serialization

using KSDFT1d

include("system_defs_01.jl")
#include("system_defs_02.jl")

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

# Read psi and Haux from a converged calculation
# Other parameters must be the same, i.e. kT
psi = deserialize("../TEMP_psi.dat")
ebands = deserialize("../TEMP_evals.dat")
Haux = diagm(0 => ebands[:,1])

# Some sanity checks
if size(psi,1) != Npoints
    error("Different: Npoints = $(Npoints), read from data $(size(psi,1))")
end

if size(psi,2) != Nstates
    error("Different: Nstates = $(Nstates), read from data $(size(psi,2))")
end

@assert size(ebands,1) == Nstates


# Evaluate total energy by calling Lfunc
E1 = calc_Lfunc_Haux!(Ham, psi, Haux)
# Check whether E1 is the same as the converged value obtained
# by SCF calculation

g = zeros(Npoints,Nstates)
Hsub = zeros(Nstates,Nstates)
g_Haux = zeros(Nstates,Nstates)
Kg_Haux = zeros(Nstates,Nstates)

# Evaluate gradients
calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)
# Check whether the g and g_Haux are close to zeros.

display(g_Haux)
display(Kg_Haux)