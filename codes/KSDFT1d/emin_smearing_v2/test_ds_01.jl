L = 10.0;
Ngrid = 51;
atoms = [-1.0, 1.0];
Z  = [1.0, 1.0];
a_soft_ee = 1.0;
a_soft_en = 1.0;
N_electrons = 2;
kT = 0.01;
Nstates = 5;
degeneracy = 2;
κ = 0.01;
maxiter = 200;
tol = 1e-10;
verbose = true;
seed = nothing;

grid = Grid1D(L=L, N=Ngrid);
V_ext = external_potential(grid, atoms, Z, a_soft_en);
precond_psi = preconditioner_psi(grid);

# initial wavefunctions
psi = randn(Ngrid, Nstates);
orthonormalize!(psi, grid.dx);

# initial diagonal η (guess energies)
ε = collect(range(-1.0, 1.0, length=Nstates));
μ = find_mu(ε, N_electrons, kT, degeneracy);
f = fermi.(ε, μ, kT, degeneracy);

# initial density and potentials
ρ = compute_density(psi, f);
V_H = hartree_potential(ρ, grid, a_soft_ee);
V_tot = V_ext + V_H;

# kinetic energy only
Tpsi = apply_kinetic(psi, grid);
Tmat = dot_grid(psi, Tpsi, grid.dx);

# full Hamiltonian (for gradients)
Hpsi = apply_H(psi, V_tot, grid);
Hmat = dot_grid(psi, Hpsi, grid.dx);

F = free_energy(psi, f, Tmat, ρ, V_H, grid, V_ext, kT, degeneracy);