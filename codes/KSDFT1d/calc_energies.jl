function calc_E_NN(atoms::Atoms1d)
    Natoms = atoms.Natoms
    Zvals = atoms.Zvals
    atpos = atoms.positions
    E_NN = 0.0
    for i in 1:Natoms, j in (i+1):Natoms
        rij = abs(atpos[i] - atpos[j])
        E_NN += Zvals[i]*Zvals[j]/rij
    end
    return E_NN
end

function calc_E_kin(Ham, psi)
    Ekin = 0.0
    Nstates = size(psi, 2)
    Kpsi = Ham.Kmat*psi
    Focc = Ham.electrons.Focc
    dx = Ham.grid.dx
    for ist in 1:Nstates
        Ekin += Focc[ist]*dot(psi[:,ist], Kpsi[:,ist])*dx
    end
    return Ekin
end