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
calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Kg, Hsub, g_Haux, Kg_Haux)

energies1 = deepcopy(Ham.energies)
println("New variables:")
for ist in 1:Nstates
    @printf("%5d %18.10f %18.10f\n", ist, Ham.electrons.Focc[ist], Ham.electrons.ebands[ist,1])
end
@printf("Ham.electrons.E_fermi = %18.10f\n", Ham.electrons.E_fermi)


# Also try to vary only psi or only Haux

# Set the directions
Kg = zeros(Float64, size(g))
prec_invK!(Ham, g, Kg)
dW = Kg # why using Kg will not work? Do we need to constrain the search direction?
constrain_search_dir!(dW, psi, hx)
dW_Haux = g_Haux
for iexp in 2:10
    #
    Δ = 10.0^(-iexp)
    println("\nStart Δ = ", Δ)
    #
    psi_new = psi + Δ*dW
    Haux_new = Haux + Δ*dW_Haux
    #
    # Prepare new psi and Haux
    Udagger = inv(sqrt(psi_new'*psi_new)) ./ sqrt(hx)
    psi_new[:,:] = psi_new*Udagger
    #Haux_new = Udagger' * Haux_new * Udagger
    Urot = transform_psi_Haux!(psi_new, Haux_new)
    #
    # Evaluate energy at new variables
    E_new = calc_Lfunc_Haux!(Ham, psi_new, Haux_new)
    println("New variables:")
    for ist in 1:Nstates
        @printf("%5d %18.10f %18.10f\n", ist, Ham.electrons.Focc[ist], Ham.electrons.ebands[ist,1])
    end
    @printf("Ham.electrons.E_fermi = %18.10f\n", Ham.electrons.E_fermi)

    energies_new = deepcopy(Ham.energies)

    println("E1    = ", E1)
    println("E_new = ", E_new)
    dE = abs(E_new - E1)
    println("dE      = ", dE)

    # XXX: Use different gradient from this? from Lfunc_Focc?
    println()
    dE_psi = 2*real(dot(g, Δ*dW)*hx)
    dE_Haux = real(dot(g_Haux, Δ*dW_Haux))
    println("dE_psi  = ", dE_psi)
    println("dE_Haux = ", dE_Haux)
    println("sum = ", dE_psi + dE_Haux)
    println("ratio = ", (dE_psi + dE_Haux)/dE)
end

