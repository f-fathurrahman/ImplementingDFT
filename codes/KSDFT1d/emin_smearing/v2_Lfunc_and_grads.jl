
# Update of Hamiltonian is done outside


# Haux is assumed to be already diagonalized
# and psi is already rotated accordingly
function update_from_wavefunc_Haux!(Ham, psi, Haux)
    ebands = reshape(diag(Haux), Ham.electrons.Nstates, 1)
    update_from_ebands!(Ham, ebands)
    update_from_wavefunc!(Ham, psi)
    return
end


function v2_calc_Lfunc_Haux!(
    Ham::Hamiltonian1d,
    psi::Matrix{Float64}, # (Nbasis,Nstates)
    Haux::Matrix{Float64} # (Nstates,Nstates)
)

    hx = Ham.grid.hx    
    Npoints = Ham.grid.Npoints
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    Vtot = Ham.potentials.Total
    rhoe = Ham.rhoe
    mTS = Ham.energies.mTS

    # Calculate ebands first
    # Because Haux is assumed to be already diagonalized, we simply
    # take the diagonal elements
    ebands = reshape(diag(Haux), Nstates, 1)
    # This is for Nspin = 1

    # Evaluate total energy components
    Ekin = calc_E_kin(Ham, psi)
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

    #println(Ham.energies)

    return Etot
end



function v2_calc_grad_Lfunc_Haux!(
    Ham::Hamiltonian1d,
    psi, # (Nbasis,Nstates)
    Haux, # (Nstates,Nstates)
    g,
    Hsub,
    g_Haux,
    Kg_Haux
)
    fill!(g, 0.0)
    fill!(Hsub, 0.0)
    fill!(g_Haux, 0.0)
    fill!(Kg_Haux, 0.0)
    # Evaluate the gradient for psi
    calc_grad!(Ham, psi, g, Hsub)
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
    return
end