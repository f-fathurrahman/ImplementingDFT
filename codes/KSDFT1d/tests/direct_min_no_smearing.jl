# Evaluate gradient at psi
# Ham matrix is not updated
function calc_grad( Ham, psi; update=false )

    if update
        update_from_wavefunc!(Ham, psi)
    end
    
    ispin = 1
    Focc = Ham.electrons.Focc
    hx = Ham.grid.hx

    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    #
    grad = zeros(Float64, Nbasis, Nstates)

    Hpsi = Ham*psi
    for i in 1:Nstates
        @views grad[:,i] .= Hpsi[:,i]
        for j in 1:Nstates
            @views λ = dot(psi[:,j], Hpsi[:,i]) * hx
            @views grad[:,i] .= grad[:,i] .-  λ * psi[:,j]
        end
        grad[:,i] .= Focc[i,ispin]*grad[:,i]
    end

    # the usual case of constant occupation numbers
    if all(Focc[:,ispin] .== 2.0) == true || all(Focc[:,ispin] .== 1.0) == true
        # immediate return
        return grad
    end

    #F = Matrix(Diagonal(Focc[:,ispin]))
    #ℍ = psi' * Hpsi * hx
    #HFH = ℍ*F - F*ℍ
    #ℚ = 0.5*HFH
    #@views grad[:,:] .+= psi*ℚ # additional contributions
    
    return grad
end

function update_from_wavefunc!(Ham, psi)
    Npoints = Ham.grid.Npoints
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    Vtot = Ham.potentials.Total
    rhoe = Ham.rhoe

    # Electron density
    calc_rhoe!(Ham, psi, rhoe)
    #println("integ rhoe = ", sum(rhoe)*hx)

    # Update the potentials
    # Note that Vxc, Vhartree, and Vtot refers to Ham.potentials
    ρ = reshape(rhoe, Npoints) # FIXME: need to do this is Poisson_solve_sum!
    Poisson_solve_sum!(Ham.grid, ρ, Vhartree)
    Vxc[:] = calc_Vxc_1d(Ham.xc_calc, rhoe)
    @views Vtot[:,1] = Vion[:] + Vhartree[:] + Vxc[:,1]

    return
end

# Ham will be modified (the potential terms)
function calc_KohnSham_Etotal!(Ham, psi)
    
    update_from_wavefunc!(Ham, psi)

    hx = Ham.grid.hx    
    Npoints = Ham.grid.Npoints
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    Vtot = Ham.potentials.Total
    rhoe = Ham.rhoe

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
    Etot = Ekin + Ehartree + Eion + Exc + Ham.energies.NN

    return Etot
end