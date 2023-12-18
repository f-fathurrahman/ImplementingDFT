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
E1 = calc_Lfunc_Haux!(Ham, psi, Haux)
calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)
println("Original variables:")
println("Ham.electrons.ebands  = ", Ham.electrons.ebands)
println("Ham.electrons.Focc    = ", Ham.electrons.Focc)
println("Ham.electrons.E_fermi = ", Ham.electrons.E_fermi)


Δ = 1e-8
Δ_Haux = 0.0
dW = randn(Float64, size(psi))
dW_Haux = randn(Float64, size(Haux))

psi_new = psi + Δ*dW
Haux_new = Haux + Δ_Haux*dW_Haux

# Prepare new psi and Haux
Udagger = inv(sqrt(psi_new'*psi_new)) ./ sqrt(hx)
psi_new[:,:] = psi_new*Udagger
Haux_new = Udagger' * Haux_new * Udagger
Urot = transform_psi_Haux!(psi_new, Haux_new)
E_new = calc_Lfunc_Haux!(Ham, psi_new, Haux_new)
println("New variables:")
println("Ham.electrons.ebands  = ", Ham.electrons.ebands)
println("Ham.electrons.Focc    = ", Ham.electrons.Focc)
println("Ham.electrons.E_fermi = ", Ham.electrons.E_fermi)

energies_new = deepcopy(Ham.energies)

println("E1    = ", E1)
println("E_new = ", E_new)
dE = abs(E_new - E1)
println("dE      = ", dE)

# XXX: Use different gradient from this? from Lfunc_Focc?
println()
dE_psi = abs( 2*real(dot(g, Δ*dW)*hx) )
dE_Haux = abs( real(dot(g_Haux, Δ_Haux*dW_Haux)) )
println("dE_psi  = ", dE_psi)
println("dE_Haux = ", dE_Haux)
println("sum = ", dE_psi + dE_Haux)
println("ratio = ", (dE_psi + dE_Haux)/dE)

println()
dE_psi = abs( 2*real(dot(g, Δ*dW*inv(Urot))*hx) )
dE_Haux = abs( real(dot(g_Haux, Δ_Haux*dW_Haux)) )
println("dE_psi  = ", dE_psi)
println("dE_Haux = ", dE_Haux)
println("sum = ", dE_psi + dE_Haux)
println("ratio = ", (dE_psi + dE_Haux)/dE)


println()
dE_psi = abs( 2*real(dot(g, psi-psi_new)*hx) )
dE_Haux = abs( real(dot(g_Haux, Haux-Haux_new)) )
println("dE_psi  = ", dE_psi)
println("dE_Haux = ", dE_Haux)
println("sum = ", dE_psi + dE_Haux)
println("ratio = ", (dE_psi + dE_Haux)/dE)
