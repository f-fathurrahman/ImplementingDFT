function calc_E_nn(atoms::Atoms1d)
    @assert atoms.pbs == false
    Natoms = atoms.Natoms
    Enn = 0.0
    Zvals = atoms.Zvals
    atpos = atoms.positions
    for ia in 1:Natoms
        for ja in 2:Natoms
            r = abs(atpos[ia] - atpos[ja]) # in 1d
            Enn += Zvals[ia]*Zvals[ja]/r
        end
    end
    Enn = 0.5*Enn # correct for double counting
    return Enn
end

function calc_E_kin(Ham, psi)
    Ekin = 0.0
    Nstates = size(psi, 2)
    Kpsi = Ham.Kmat*psi
    Focc = Ham.electrons.Focc
    hx = Ham.grid.hx
    for ist in 1:Nstates
        Ekin += Focc[ist]*dot(psi[:,ist], Kpsi[:,ist])*hx
    end
    return Ekin
end