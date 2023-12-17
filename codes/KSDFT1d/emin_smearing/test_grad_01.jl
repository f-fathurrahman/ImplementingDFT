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
include("v2_Lfunc_and_grads.jl")

Ham = init_Hamiltonian()

hx = Ham.grid.hx
Npoints = Ham.grid.Npoints
Nstates = Ham.electrons.Nstates

# Random wavefunc
#Random.seed!(1234)
psi = generate_random_wavefunc(Ham)
psi_orig = deepcopy(psi)

Haux = diagm( 0 => sort(randn(Nstates)) )
Haux_orig = deepcopy(Haux)

# psi and Haux are ready

# Allocate memory
g = zeros(Float64, size(psi))
g_Haux = zeros(Float64, size(Haux))
Kg_Haux = zeros(Float64, size(Haux))
Hsub = zeros(Float64, size(Haux))

# Evaluate total energy by calling Lfunc
update_from_wavefunc_Haux!(Ham, psi, Haux)
E1 = v2_calc_Lfunc_Haux!(Ham, psi, Haux)
energies1 = deepcopy(Ham.energies)
v2_calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)

Δ = 1e-5
dW = randn(Float64, size(psi))
dW_Haux = randn(Float64, size(Haux))

psi_new = psi + Δ*dW
Haux_new = Haux + Δ*dW_Haux

# Orthonormalize (involves rotation)
Udagger = inv(sqrt(psi_new'*psi_new)) ./ sqrt(hx)
psi_new[:,:] = psi_new*Udagger
Haux_new = Udagger' * Haux_new * Udagger
Urot = transform_psi_Haux!(psi_new, Haux_new)


update_from_wavefunc_Haux!(Ham, psi_new, Haux_new)
E_new = v2_calc_Lfunc_Haux!(Ham, psi_new, Haux_new)

energies_new = deepcopy(Ham.energies)

println("E1    = ", E1)
println("E_new = ", E_new)

dE = abs(E_new - E1)
println("dE      = ", dE)

dE_psi = abs( 2*real(dot(g, Δ*dW)*hx) )
dE_Haux = abs( real(dot(g_Haux, Δ*dW_Haux)) )
println("dE_psi  = ", dE_psi)
println("dE_Haux = ", dE_Haux)
println("sum = ", dE_psi + dE_Haux)

println("ratio = ", (dE_psi + dE_Haux)/dE)