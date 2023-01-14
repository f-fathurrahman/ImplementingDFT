function update_from_Haux!(Ham, Haux)
    kT = 0.1*eV2Ha # XXX FIXED
    Ham.electrons.kT = kT

    # Calculate Kohn-Sham eigenvalues and occupation numbers
    Focc = Ham.electrons.Focc
    ebands = Ham.electrons.ebands
    
    println("Debug Haux: ", Haux[14,14])

    ebands[:,1], Urot = eigen(Hermitian(Haux)) # need to symmetrize?
    Nelectrons = Ham.electrons.Nelectrons
    Ham.electrons.E_fermi, Ham.energies.mTS = update_Focc!(
        Focc, smear_fermi, smear_fermi_entropy,
        ebands, Float64(Nelectrons), kT
    )

    println("---------------------")
    println("In update_from_Haux: ")
    println("---------------------")
    println()
    println("ebands = ", ebands)
    println("E_fermi = ", Ham.electrons.E_fermi)
    println("mTS = ", Ham.energies.mTS)
    println()

    return Urot
end


function update_from_wavefunc!(Ham, psi)
    Npoints = Ham.grid.Npoints
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    Vtot = Ham.potentials.Total
    rhoe = Ham.rhoe

    println("Enter update_from_wavefunc!:")

    # Electron density
    calc_rhoe!(Ham, psi, rhoe)
    println("integ rhoe = ", sum(rhoe)*hx)

    # Update the potentials
    # Note that Vxc, Vhartree, and Vtot refers to Ham.potentials
    ρ = reshape(rhoe, Npoints) # FIXME: need to do this is Poisson_solve_sum!
    Poisson_solve_sum!(Ham.grid, ρ, Vhartree)
    Vxc[:] = calc_Vxc_1d(Ham.xc_calc, rhoe)
    @views Vtot[:,1] = Vion[:] + Vhartree[:] + Vxc[:,1]

    return
end

# Ham will be modified (the potential terms)
function calc_KohnSham_Etotal!(
    Ham,
    psi::Matrix{Float64},
    Haux::Matrix{Float64}
)

    Urot = update_from_Haux!(Ham, Haux)
    update_from_wavefunc!(Ham, psi*Urot)

    hx = Ham.grid.hx    
    Npoints = Ham.grid.Npoints
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    Vtot = Ham.potentials.Total
    rhoe = Ham.rhoe
    mTS = Ham.energies.mTS

    # Evaluate total energy components
    Ekin = calc_E_kin(Ham, psi*Urot) # XXX: use psi*Urot ?
    Ham.energies.Kinetic = Ekin
    
    Ehartree = 0.5*dot(rhoe[:,1], Vhartree)*hx
    Ham.energies.Hartree = Ehartree

    Eion = dot(rhoe, Vion)*hx
    Ham.energies.Ion = Eion

    epsxc = calc_epsxc_1d(Ham.xc_calc, rhoe[:,1])
    Exc = dot(rhoe, epsxc)*hx
    Ham.energies.XC = Exc

    # The total energy (also include nuclei-nuclei or ion-ion energy)
    Etot = Ekin + Ehartree + Eion + Exc + Ham.energies.NN + mTS

    println(Ham.energies)

    return Etot
end