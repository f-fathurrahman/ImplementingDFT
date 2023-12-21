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
include("linemin_quad.jl")
include("linemin_armijo.jl")
include("linemin_grad.jl")

Ham = init_Hamiltonian()

hx = Ham.grid.hx
Npoints = Ham.grid.Npoints
Nstates = Ham.electrons.Nstates

# Random wavefunc
iseed = abs(rand(Int64))
println("iseed = ", iseed)
Random.seed!(iseed)
psi = generate_random_wavefunc(Ham)
psi_orig = deepcopy(psi)

Haux = diagm( 0 => sort(randn(Nstates)) )
Haux_orig = deepcopy(Haux)

# psi and Haux are ready

# Allocate memory
g = zeros(Float64, size(psi))
Kg = zeros(Float64, size(psi))
g_Haux = zeros(Float64, size(Haux))
Kg_Haux = zeros(Float64, size(Haux))
Hsub = zeros(Float64, size(Haux))

# Evaluate total energy by calling Lfunc
E1 = calc_Lfunc_Haux!(Ham, psi, Haux)
calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)

prec_invK!(Ham, g, Kg) # Precondition
#calc_grad_no_Focc!(Ham, psi, Kg)
#prec_invK!(Ham, Kg)

# Set direction
d = -Kg
d_Haux = -Kg_Haux
constrain_search_dir!(d, psi, hx)

α, is_linmin_success = linemin_quad(Ham, psi, Haux, g, g_Haux, d, d_Haux, E1)

#α = linemin_grad(Ham, psi, Haux, g, g_Haux, d, d_Haux, E1)
#is_linmin_success = true # force to true

println("α = ", α)
ΔEdir = 2*dot(g, α*d)*hx + dot(g_Haux, α*d_Haux)
println("Expected ΔE = ", ΔEdir)

# Do the step
psi_new = psi + α*d
Haux_new = Haux + α*d_Haux
Udagger, Urot = prepare_psi_Haux!(psi_new, Haux_new, hx)
# Evaluate at the new variables
E_new = calc_Lfunc_Haux!(Ham, psi_new, Haux_new)
calc_grad_Lfunc_Haux!(Ham, psi_new, Haux_new, g, Hsub, g_Haux, Kg_Haux)

println("E_new = ", E_new)
if E_new > E1
    println("!!! WARNING: E_new is higher than E1")
end
println("E_new - E1 = ", E_new - E1)


