# Various functions to update Hamiltonian


# Input: ebands
# Modifies: ebands (copy), Focc, E_fermi, mTS
# Also set kT (hardcoded)
function update_from_ebands!(Ham, ebands)
    
    println("\nENTER update_from_ebands!")

    kT = 0.1*eV2Ha # XXX FIXED
    Ham.electrons.kT = kT

    # Calculate Kohn-Sham eigenvalues and occupation numbers
    Focc = Ham.electrons.Focc
    Ham.electrons.ebands[:] .= ebands[:]

    Nelectrons = Ham.electrons.Nelectrons
    Ham.electrons.E_fermi, Ham.energies.mTS = update_Focc!(
        Focc, smear_fermi, smear_fermi_entropy,
        ebands, Float64(Nelectrons), kT
    )

    println("-----------------------")
    println("In update_from_ebands!:")
    println("-----------------------")
    println()
    println("E_fermi = ", Ham.electrons.E_fermi)
    println("mTS = ", Ham.energies.mTS)
    for ist in 1:Nstates
        @printf("%3d %18.10f %18.10f\n", ist, Focc[ist,1], ebands[ist,1])
    end
    println("sum Focc: ", sum(Focc))
    println()

    println("EXIT update_from_ebands!\n")

    return
end


# Input: psi
# Modifies: Rhoe, potentials, energies
function update_from_wavefunc!(Ham, psi)
    Npoints = Ham.grid.Npoints
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    Vtot = Ham.potentials.Total
    rhoe = Ham.rhoe

    println("\nENTER update_from_wavefunc!")

    # Electron density
    calc_rhoe!(Ham, psi, rhoe)
    println("integ rhoe = ", sum(rhoe)*hx)

    # Update the potentials
    # Note that Vxc, Vhartree, and Vtot refers to Ham.potentials
    ρ = reshape(rhoe, Npoints) # FIXME: need to do this is Poisson_solve_sum!
    Poisson_solve_sum!(Ham.grid, ρ, Vhartree)
    Vxc[:] = calc_Vxc_1d(Ham.xc_calc, rhoe)
    @views Vtot[:,1] = Vion[:] + Vhartree[:] + Vxc[:,1]

    println("EXIT update_from_wavefunc!\n")

    return
end