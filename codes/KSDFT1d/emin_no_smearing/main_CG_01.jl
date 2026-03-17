# run setup_first

# can copy-pasted to REPL
Ham = init_Hamiltonian();
psi = generate_random_wavefunc(Ham);
solve_emin_CG_v01!(Ham, psi=psi);

