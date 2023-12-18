function create_atoms()
    Natoms = 2
    σ = ones(Float64, Natoms)*(1.0)
    masses = ones(Float64, Natoms)*42000.0
    Zvals = ones(Float64, Natoms)*2
    L = 10.0
    atpos = zeros(Float64, Natoms)
    atpos[1] = -1.0
    atpos[2] =  1.0
    return Atoms1d( atpos, Zvals, σ, masses, L )
end

function pot_gaussian( x, x0 )
    return -exp(-(x - x0)^2)/sqrt(1.0/π)^3
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
    Ham = Hamiltonian1d(atoms, 51)
    init_Vions!(Ham)
    return Ham
end