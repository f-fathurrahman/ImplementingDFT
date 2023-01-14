include("test_direct_min_02.jl")
include("gradients_psi_Haux.jl")

Ham = init_Hamiltonian()
Ham.energies.NN = calc_E_NN(Ham.atoms)

hx = Ham.grid.hx
Npoints = Ham.grid.Npoints
Nelectrons = Ham.electrons.Nelectrons
Nstates = Ham.electrons.Nstates
Focc = Ham.electrons.Focc

Random.seed!(1234)
psi = generate_random_wavefunc(Ham)

Haux = psi' * (Ham*psi) * hx # Hsub, subspace Hamiltonian
#Haux = rand(Nstates,Nstates)
#Haux = 0.5*(Haux + Haux')

E1 = calc_KohnSham_Etotal!(Ham, psi, Haux)

Haux_orig = copy(Haux)
Δ = 0.01

idx = CartesianIndex(14,13)
Haux[idx] = Haux_orig[idx] + Δ
Ep = calc_KohnSham_Etotal!(Ham, psi, Haux)

Haux[idx] = Haux_orig[idx] - Δ
Em = calc_KohnSham_Etotal!(Ham, psi, Haux)

g_Haux_fd = (Ep - Em)/(2*Δ)

g = zeros(Npoints,Nstates)
Hsub = zeros(Nstates,Nstates)
calc_grad!(Ham, psi, g, Hsub)

g_Haux = zeros(Nstates,Nstates)
Kg_Haux = zeros(Nstates,Nstates)
calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)