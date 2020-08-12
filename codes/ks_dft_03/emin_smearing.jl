function compute!(
    Ham::Hamiltonian, evars::ElecVars,
    g::ElecGradient, Kg::ElecGradient, kT::Float64,
    subrot::SubspaceRotations
)

    Etot = calc_energies_grad!( Ham, evars, g, Kg, kT )
    
    rotPrevCinv = subrot.prevCinv
    rotPrev = subrot.prev
    
    g.psi = g.psi * rotPrevCinv
    Kg.psi = Kg.psi * rotPrevCinv
    
    g.Haux = rotPrev * g.Haux * rotPrev'
    Kg.Haux = rotPrev * Kg.Haux * rotPrev'

    # No caching is done (for SubspaceRotationAdjutst)
    return Etot
end


function calc_energies_grad!(
    Ham::Hamiltonian, evars::ElecVars,
    g::ElecGradient, Kg::ElecGradient, kT::Float64
)

    println("\nCalculating energy and gradients\n")

    Nstates = Ham.electrons.Nstates

    # Using electrons.eorbs as energy eigevalues
    E_f, Ham.energies.mTS =
    update_Focc!( Ham.electrons.Focc, smear_fermi, smear_fermi_entropy,
                  Ham.electrons.eorbs, Float64(Ham.electrons.Nelectrons), kT )

    @printf("Updating Focc: Fermi energy = %18.10f\n\n", E_f)

    Rhoe = zeros(Float64, Ham.grid.Npoints)
    calc_rhoe!( Ham, evars.psi, Rhoe )
    update!( Ham, Rhoe )
    
    # Calculate total energy
    calc_energies!( Ham, evars.psi )

    #
    # Gradient for psi
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
    
    dmuNum = w*sum(fprimeNum)
    dmuDen = w*sum(fprime)
    dmuContrib = dmuNum/dmuDen
    
    gradF0 = evars.Hsub - diagm( 0 => Ham.electrons.eorbs )
        
    gradF = copy(gradF0)
    for ist in 1:Nstates
        gradF[ist,ist] = gradF0[ist,ist] - dmuContrib
    end
    g_tmp = grad_smear( smear_fermi, smear_fermi_prime, Ham.electrons.eorbs, E_f, kT, gradF )
    
    g.Haux = w * 0.5 * (g_tmp' + g_tmp) # 
    Kg.Haux = -copy(gradF0) #-0.1*copy(gradF0), precondition?

    println("Finished calculating energy and gradients\n")

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
    Hsub[:,:] = ψ' * Hψ * Ham.grid.dVol
    Hψ = Hψ - ψ*Hsub
    #
    for ist in 1:Nstates
        @views g[:,ist] = Focc[ist] * Hψ[:,ist]
    end
    return
end



function do_step!(
    Ham::Hamiltonian,
    α::Float64, α_Haux::Float64, evars::ElecVars, d::ElecGradient,
    subrot::SubspaceRotations
)
    dVol = Ham.grid.dVol
    Nstates = size(evars.psi,2)

    rotPrev = subrot.prev
    rotPrevC = subrot.prevC
    rotPrevCinv = subrot.prevCinv

    evars.psi = evars.psi + α*d.psi*rotPrevC

    # Haux fillings:
    Haux = diagm( 0 => Ham.electrons.eorbs )
    Haux = Haux + α_Haux*( rotPrev' * d.Haux * rotPrev )

    Ham.electrons.eorbs, rot = eigen(Hermitian(Haux)) # need to symmetrize?
 
    #rotC = rot
    #eVars.orthonormalize(q, &rotC);
    Udagger = inv( sqrt( evars.psi' * evars.psi * dVol ) )
    rotC = Udagger*rot
    evars.psi = evars.psi*rotC

    # Accumulate rotation
    subrot.prev = rotPrev * rot
    subrot.prevC = rotPrevC * rotC
    subrot.prevCinv = inv(rotC) * rotPrevCinv
    return 
end


function constrain_search_dir!( d::ElecGradient, evars::ElecVars, dVol )
    d.psi = d.psi - evars.psi * ( evars.psi' * d.psi * dVol )
    return
end