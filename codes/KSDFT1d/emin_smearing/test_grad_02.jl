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
println("Before calculating energy:")
println("ebands = ", Ham.electrons.ebands)
println("ebands = ", Ham.electrons.Focc)
E0 = calc_Lfunc_Haux!(Ham, psi, Haux)
println("After calculating energy:")
println("ebands = ", Ham.electrons.ebands)
println("ebands = ", Ham.electrons.Focc)
calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)
println("After calculating gradients:")
println("ebands = ", Ham.electrons.ebands)
println("ebands = ", Ham.electrons.Focc)


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
E_plus = calc_Lfunc_Haux!(Ham, psi_new, Haux_new)
println("After calculating energy:")
println("ebands = ", Ham.electrons.ebands)
println("ebands = ", Ham.electrons.Focc)


#=
dW[1,1] = -Δ
psi_new = psi + dW
Haux_new = Haux + dW_Haux
# Orthonormalize (involves rotation)
Udagger = inv(sqrt(psi_new'*psi_new)) ./ sqrt(hx)
psi_new[:,:] = psi_new*Udagger
Haux_new = Udagger' * Haux_new * Udagger
Urot = transform_psi_Haux!(psi_new, Haux_new)
E_minus = calc_Lfunc_Haux!(Ham, psi_new, Haux_new)
=#
