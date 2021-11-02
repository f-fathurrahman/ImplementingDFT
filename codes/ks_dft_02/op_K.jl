# Kinetic operator
function op_K( Ham::Hamiltonian, psi )
    return -0.5*Ham.∇2 * psi
end