function create_atoms()
    Natoms = 1
    σ = ones(Float64, Natoms)*(1.0)
    masses = ones(Float64, Natoms)*42000.0
    Zvals = ones(Float64, Natoms)*3.0
    L = 10.0
    atpos = zeros(Float64, Natoms)
    atpos[1] = 0.0
    return Atoms1d( atpos, Zvals, σ, masses, L )
end

function pot_gaussian( x, x0 )
    return -5.0*exp(-4.5*(x-x0)^2)
end

function init_Vions!(Ham)
    atpos = Ham.atoms.positions
    Natoms = Ham.atoms.Natoms
    for ia in 1:Natoms
        Ham.potentials.Ions[:] += pot_gaussian.(Ham.grid.x, atpos[ia])
    end
    return
end

function init_Hamiltonian()
    atoms = create_atoms()
    Ham = Hamiltonian1d(atoms, 51, Nstates_extra=5)
    init_Vions!(Ham)
    Ham.energies.NN = calc_E_NN(Ham.atoms) # also calculate E_NN
    return Ham
end