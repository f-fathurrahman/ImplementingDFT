function create_atoms_03()
    Natoms = 2
    σ = ones(Float64, Natoms)*(1.0)
    masses = ones(Float64, Natoms)*42000.0
    Zvals = ones(Float64, Natoms)*5.0
    L = 10.0
    atpos = zeros(Float64, Natoms)
    atpos[1] = -1.5
    atpos[2] =  1.5
    return Atoms1d( atpos, Zvals, σ, masses, L )
end

function pot_gaussian_03( x, x0 )
    return -25.0*exp(-4.5*(x-x0)^2)
end
# too negative amplitude might gave NaN to g_Haux

function init_Vions_03!(Ham)
    atpos = Ham.atoms.positions
    Natoms = Ham.atoms.Natoms
    for ia in 1:Natoms
        Ham.potentials.Ions[:] += pot_gaussian_03.(Ham.grid.x, atpos[ia])
    end
    return
end

function init_Hamiltonian_03()
    atoms = create_atoms_03()
    Ham = Hamiltonian1d(atoms, 61, Nstates_extra=6)
    init_Vions_03!(Ham)
    Ham.electrons.kT = 0.1*eV2Ha
    Ham.energies.NN = calc_E_NN(Ham.atoms) # also calculate E_NN
    return Ham
end