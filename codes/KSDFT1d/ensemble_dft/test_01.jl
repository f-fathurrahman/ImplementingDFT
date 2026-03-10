function create_atoms_01()
    Natoms = 2
    σ = ones(Float64, Natoms)*(1.0)
    masses = ones(Float64, Natoms)
    #Zvals = ones(Float64, Natoms)
    Zvals = [1.0, 1.0]
    L = 10.0  # from -L/2 to +L/2
    atpos = zeros(Float64, Natoms)
    atpos[1] = -2.0
    atpos[2] =  2.0
    return Atoms1d( atpos, Zvals, σ, masses, L )
end

function pot_soft_coulomb( x, x0; Z=1.0, a=1.0 )
    return -Z/sqrt((x - x0)^2 + a^2)
end

function init_Vions_01!(Ham)
    atpos = Ham.atoms.positions
    Natoms = Ham.atoms.Natoms
    Zvals = Ham.atoms.Zvals
    for ia in 1:Natoms
        Ham.potentials.Ions[:] += pot_soft_coulomb.(Ham.grid.x, atpos[ia], Z=Zvals[ia])
    end
    return
end

function init_Hamiltonian_01()
    atoms = create_atoms_01()
    Ham = Hamiltonian1d(atoms, 51, Nstates_extra=6)
    init_Vions_01!(Ham)
    Ham.electrons.kT = 0.01 #*eV2Ha
    Ham.energies.NN = calc_E_NN(Ham.atoms) # also calculate E_NN
    return Ham
end


