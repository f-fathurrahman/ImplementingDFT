# Kinetic operator
function op_K( Ham::Hamiltonian, psi )
    return -0.5*Ham.Laplacian * psi
end