push!(LOAD_PATH, "../")

import Random
using Printf
using LinearAlgebra
using Serialization

using KSDFT1d

#include("system_defs_01.jl")
#include("system_defs_02.jl")
include("system_defs_03.jl")

include("Lfunc.jl")
include("../utilities.jl")
include("gradients_psi_Haux.jl")

Ham = init_Hamiltonian()

hx = Ham.grid.hx
Npoints = Ham.grid.Npoints
Nstates = Ham.electrons.Nstates

# Random wavefunc
#Random.seed!(111) # vary this to find problematic case?
psi = generate_random_wavefunc(Ham)
Haux = diagm( 0 => sort(randn(Nstates)) )
#
d = zeros(Float64, size(psi))
d_Haux = zeros(Float64, size(Haux))
# Allocate memory
g = zeros(Float64, size(psi))
Kg = zeros(Float64, size(psi))
#
g_Haux = zeros(Float64, size(Haux))
Kg_Haux = zeros(Float64, size(Haux))
Hsub = zeros(Float64, size(Haux))

E1 = calc_Lfunc_Haux!(Ham, psi, Haux)
calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)
#prec_invK!(Ham, g, Kg)
calc_grad_no_Focc!(Ham, psi, Kg)
prec_invK!(Ham, Kg)

@printf("%18.10f %18.10f\n", 0.0, E1)

α_start = 1e-10
mult_factor = 5
Nsteps = 10
α = α_start
for i in 1:Nsteps
    α *= mult_factor
    d[:,:] = -g
    d_Haux[:,:] = -Kg_Haux
    #
    psi_new = psi + α*d
    Haux_new = Haux + α*d_Haux
    prepare_psi_Haux!(psi_new, Haux_new, hx)
    E_new = calc_Lfunc_Haux!(Ham, psi_new, Haux_new)
    dg = 2*dot(g, α*d)*hx + dot(g_Haux, α*d_Haux)
    dE = E_new - E1
    ratio = dE/dg
    @printf("%18.10f %18.10f %18.10f %18.10f %18.10f\n", α, E_new, dE, dg, ratio)
end

