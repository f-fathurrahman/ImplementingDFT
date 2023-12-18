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
update_from_wavefunc!(Ham, psi) # using default Focc

# Calculate band energies from subspace Hamiltonian obtained from current psi
Hsub = psi' * (Ham * psi) * hx
ebands = reshape(eigvals(Hsub), Nstates, 1)
Haux = diagm( 0 => ebands[:,1] )
Haux_orig = deepcopy(Haux)

# psi and Haux are ready

# Allocate memory
g = zeros(Float64, size(psi))
g_Haux = zeros(Float64, size(Haux))
Kg_Haux = zeros(Float64, size(Haux))
Hsub = zeros(Float64, size(Haux))

update_from_ebands!(Ham, ebands)

#Ham.electrons.ebands[:,:] = ebands[:,:]
## Also need to set E_fermi accordingly, probaly call to update_from_ebands?
#Nstates_occ = Ham.electrons.Nstates_occ
#println("Nstates_occ = ", Nstates_occ)
#Ham.electrons.E_fermi = ebands[Nstates_occ,1]

# Evaluate total energy by calling Lfunc
println("Before calculating energy:")
ebands_orig = deepcopy(ebands)
Focc_orig = deepcopy(Ham.electrons.Focc)

println("ebands = ", Ham.electrons.ebands)
println("Focc = ", Ham.electrons.Focc)

E0 = v2_calc_Lfunc_Haux!(Ham, psi, Haux)

println("After calculating energy:")
println("ebands = ", Ham.electrons.ebands)
println("Focc = ", Ham.electrons.Focc)

v2_calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)
println("After calculating gradients:")
println("ebands = ", Ham.electrons.ebands)
println("ebands = ", Ham.electrons.Focc)

#=
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
=#


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
