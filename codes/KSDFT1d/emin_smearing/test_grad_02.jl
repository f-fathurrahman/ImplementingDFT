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
Random.seed!(1234)
psi = generate_random_wavefunc(Ham)

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
E0 = v2_calc_Lfunc_Haux!(Ham, psi, Haux)
v2_calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)


Δ = 1e-5
dW = zeros(Float64, size(psi))
dW_Haux = zeros(Float64, size(Haux))

dW[1,1] = Δ
psi_new = psi + dW
Haux_new = Haux + dW_Haux
# Orthonormalize (involves rotation)
Udagger = inv(sqrt(psi_new'*psi_new)) ./ sqrt(hx)
psi_new[:,:] = psi_new*Udagger
Haux_new = Udagger' * Haux_new * Udagger
Urot = transform_psi_Haux!(psi_new, Haux_new)
update_from_wavefunc_Haux!(Ham, psi_new, Haux_new)
E_plus = v2_calc_Lfunc_Haux!(Ham, psi_new, Haux_new)


dW[1,1] = -Δ
psi_new = psi + dW
Haux_new = Haux + dW_Haux
# Orthonormalize (involves rotation)
Udagger = inv(sqrt(psi_new'*psi_new)) ./ sqrt(hx)
psi_new[:,:] = psi_new*Udagger
Haux_new = Udagger' * Haux_new * Udagger
Urot = transform_psi_Haux!(psi_new, Haux_new)
update_from_wavefunc_Haux!(Ham, psi_new, Haux_new)
E_minus = v2_calc_Lfunc_Haux!(Ham, psi_new, Haux_new)

