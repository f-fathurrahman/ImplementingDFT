function create_atoms_01()
    Natoms = 3
    σ = ones(Float64, Natoms)*(1.0)
    masses = ones(Float64, Natoms)*42000.0
    Zvals = ones(Float64, Natoms)*5.0
    L = 10.0
    atpos = zeros(Float64, Natoms)
    atpos[1] = -1.5
    atpos[2] =  0.1
    atpos[3] =  1.5
    return Atoms1d( atpos, Zvals, σ, masses, L )
end

function pot_gaussian_01( x, x0 )
    return -25.0*exp(-4.5*(x-x0)^2)
end

function init_Vions_01!(Ham)
    atpos = Ham.atoms.positions
    Natoms = Ham.atoms.Natoms
    for ia in 1:Natoms
        Ham.potentials.Ions[:] += pot_gaussian_01.(Ham.grid.x, atpos[ia])
    end
    return
end

function init_Hamiltonian_01()
    atoms = create_atoms_01()
    Ham = Hamiltonian1d(atoms, 51, Nstates_extra=6)
    init_Vions_01!(Ham)
    Ham.electrons.kT = 0.1*eV2Ha
    Ham.energies.NN = calc_E_NN(Ham.atoms) # also calculate E_NN
    return Ham
end