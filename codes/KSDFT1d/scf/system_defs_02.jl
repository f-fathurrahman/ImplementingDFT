function init_Hamiltonian_02(; Npoints = 51)
    
    function create_atoms()
        Natoms = 3
        σ = ones(Float64, Natoms)*(1.0)
        masses = ones(Float64, Natoms)
        #Zvals = ones(Float64, Natoms)
        Zvals = [1.0, 1.0, 1.0]
        L = 10.0  # from -L/2 to +L/2
        atpos = zeros(Float64, Natoms)
        atpos[1] = -1.0
        atpos[2] =  0.0
        atpos[3] =  1.0
        return Atoms1d( atpos, Zvals, σ, masses, L )
    end

    function pot_soft_coulomb( x, x0; Z=1.0, a=1.0 )
        return -Z/sqrt((x - x0)^2 + a^2)
    end

    function init_Vions!(Ham)
        atpos = Ham.atoms.positions
        Natoms = Ham.atoms.Natoms
        Zvals = Ham.atoms.Zvals
        for ia in 1:Natoms
            Ham.potentials.Ions[:] += pot_soft_coulomb.(Ham.grid.x, atpos[ia], Z=Zvals[ia])
        end
        return
    end

    atoms = create_atoms()
    Ham = Hamiltonian1d(atoms, Npoints, use_smearing = true)
    init_Vions!(Ham)
    Ham.electrons.kT = 0.1
    Ham.energies.NN = calc_E_NN(Ham.atoms) # also calculate E_NN
    return Ham
end


