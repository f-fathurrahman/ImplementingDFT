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

Ham = init_Hamiltonian()

hx = Ham.grid.hx
Npoints = Ham.grid.Npoints

# Random wavefunc
Random.seed!(1234)
psi = generate_random_wavefunc(Ham)

update_from_wavefunc!(Ham, psi) # update the potential
    
# Prepare diagonal Haux
ebands1 = sort(randn(Nstates))
Haux = diagm(0 => ebands1)

# Evaluate total energy by calling Lfunc
E1 = calc_Lfunc_Haux!(Ham, psi, Haux)
calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)

Δ = 1e-10
dW = randn(Float64, Npoints, Nstates)

psi1 = psi + Δ*dW

# Orthonormalize (involves rotation)
Udagger = inv(sqrt(psi1'*psi1)) ./ sqrt(hx)
psi[:,:] = psi*Udagger

Haux2 = Udagger' * Haux * Udagger

# Not sure about this ???

#E1 = calc_KohnSham_Etotal!(Ham, psi1)
#dE = 2*real( sum( g0 .* (Δ*dW) ) )*hx
#println("dE = ", dE)
#println("E1 - E0 = ", E1 - E0)
#println("ratio (should be close to 1) = ", abs((E1 - E0)/dE))