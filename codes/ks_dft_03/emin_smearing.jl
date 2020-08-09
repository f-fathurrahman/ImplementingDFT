function calc_energies_grad!(
    Ham::Hamiltonian, evars::ElecVars,
    g::ElecGradient, Kg::ElecGradient, kT::Float64
)

    println("eorbs = ", Ham.electrons.eorbs)

    # Using electrons.eorbs as energy eigevalues
    E_f, Ham.energies.mTS =
    update_Focc!( Ham.electrons.Focc, smear_fermi, smear_fermi_entropy,
                  Ham.electrons.eorbs, Float64(Ham.electrons.Nelectrons), kT )
    @printf("Fermi energy = %18.10f\n", E_f)

    Rhoe = zeros(Float64, Ham.grid.Npoints)
    calc_rhoe!( Ham, evars.psi, Rhoe )
    update!( Ham, Rhoe )
    
    # Calculate total energy
    calc_energies!( Ham, evars.psi )

    Nstates = Ham.electrons.Nstates

    #
    # Gradient for psiks
    #
    calc_grad!( Ham, evars.psi, g.psi, evars.Hsub )
    # apply preconditioner
    Kg.psi[:,:] = g.psi[:,:]
    for ist in 1:Nstates
        @views ldiv!(Ham.precKin, Kg.psi[:,ist])
    end

    #
    # Gradient for Haux
    #
    fprime = zeros(Float64,Nstates)
    fprimeNum = zeros(Float64,Nstates)
    dmuNum = 0.0  # should depend on spin
    dmuDen = 0.0  # should depend on spin

    w = 2.0 # Doubly occupied orbitals
    for ist in 1:Nstates
        fprime[ist] = smear_fermi_prime( Ham.electrons.eorbs[ist], E_f, kT )
        fprimeNum[ist] = fprime[ist] * ( evars.Hsub[ist,ist] - Ham.electrons.eorbs[ist] )
    end
    
    dmuNum = dmuNum + w*sum(fprimeNum)
    dmuDen = dmuDen + w*sum(fprime)
    dmuContrib = dmuNum/dmuDen
    
    gradF0 = evars.Hsub - diagm( 0 => Ham.electrons.eorbs )
    #gradF0[:] = diagm( 0 => ( diag(evars.Hsub[i]) - evars.Haux_eigs[:,i] ) )
        
    gradF = copy(gradF0)
    for ist in 1:Nstates
        gradF[ist,ist] = gradF0[ist,ist] - dmuContrib
    end
    g_tmp = grad_smear( smear_fermi, smear_fermi_prime, Ham.electrons.eorbs, E_f, kT, gradF )
    g.Haux = w * 0.5 * (g_tmp' + g_tmp) # 
    Kg.Haux = -copy(gradF0) #-0.1*copy(gradF0), precondition?

    # Return total free energy
    return sum( Ham.energies )
end



# Calculate gradient and subspace Hamiltonian for a given wave function.
# Hamiltonian is not updated.
function calc_grad!(
    Ham::Hamiltonian,
    ψ::Matrix{Float64},
    g::Matrix{Float64},
    Hsub::Matrix{Float64},
)
    #
    Nstates = size(ψ,2)
    Focc = Ham.electrons.Focc
    #
    Hψ = op_H( Ham, ψ )
    Hsub[:] = ψ' * Hψ * Ham.grid.dVol
    Hψ = Hψ - ψ*Hsub
    #
    for ist in 1:Nstates
        g[:,ist] = Focc[ist] * Hψ[:,ist]
    end
    return
end