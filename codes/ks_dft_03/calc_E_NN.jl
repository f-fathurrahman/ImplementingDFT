function calc_E_NN( atoms::Atoms, pspots::Vector{PsPot_GTH} )

    r = atoms.positions
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    Zvals = get_Zvals(pspots) # a safer way

    if atoms.pbc == (true,true,true)
        return calc_E_NN( atoms.LatVecs, atoms, Zvals )
    end

    E_NN = 0.0
    for ia in 1:Natoms
        for ja in ia+1:Natoms
            #
            isp = atm2species[ia]
            jsp = atm2species[ja]
            #
            dx = r[1,ja] - r[1,ia]
            dy = r[2,ja] - r[2,ia]
            dz = r[3,ja] - r[3,ia]
            #
            r_ij = sqrt(dx^2 + dy^2 + dz^2)
            #
            E_NN = E_NN + Zvals[isp]*Zvals[jsp]/r_ij
        end
    end

    return E_NN
end