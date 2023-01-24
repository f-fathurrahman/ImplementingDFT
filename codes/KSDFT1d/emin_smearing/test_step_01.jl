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

Random.seed!(1235)
psi = generate_random_wavefunc(Ham)

# Prepare Haux
Haux = psi' * (Ham*psi) * hx # Hsub, subspace Hamiltonian
# Using diagonal Haux
#ebands1 = sort(randn(Nstates))
#Haux = diagm(0 => ebands1)

Urot = transform_psi_Haux!(psi, Haux)


# Evaluate total energy by calling Lfunc
E1 = calc_Lfunc_Haux!(Ham, psi, Haux)

g = zeros(Npoints,Nstates)
Kg = zeros(Npoints,Nstates)

Hsub = zeros(Nstates,Nstates)
g_Haux = zeros(Nstates,Nstates)
Kg_Haux = zeros(Nstates,Nstates)

d = zeros(Npoints,Nstates)
d_Haux = zeros(Nstates,Nstates)

# Evaluate gradients
calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)

# Precondition
prec_invK!(Ham, g, Kg)

display(g_Haux)
display(Kg_Haux)

d[:,:] = -Kg[:,:]
d_Haux[:,:] = -Kg_Haux[:,:]

α = 3e-5 # 1e-5
α_Haux = 1e-1 # 1e-5


psi[:,:] = psi[:,:] + α*d[:,:]
Haux[:,:] = Haux[:,:] + α_Haux*d_Haux[:,:]

# Orthonormalize
ortho_sqrt!(psi)
psi[:,:] = psi[:,:]*( 1.0/sqrt(hx) )


Urot2 = transform_psi_Haux!(psi, Haux)


E2 = calc_Lfunc_Haux!(Ham, psi, Haux)
calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)

println("E1 = ", E1)
println("E2 = ", E2)

if E1 < E2
    println("WARNING: Energy does not decrease")
end
