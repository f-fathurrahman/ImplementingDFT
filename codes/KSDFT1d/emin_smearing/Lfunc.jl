mutable struct RotationsCache
    rotPrev::Vector{Matrix{Float64}}
    rotPrevC::Vector{Matrix{Float64}}
    rotPrevCinv::Vector{Matrix{Float64}}
    Urot::Vector{Matrix{Float64}}
    UrotC::Vector{Matrix{Float64}}
end

function RotationsCache(Nspin, Nstates)
    rotPrev = Vector{Matrix{Float64}}(undef, Nspin)
    rotPrevC = Vector{Matrix{Float64}}(undef, Nspin)
    rotPrevCinv = Vector{Matrix{Float64}}(undef, Nspin)
    Urot = Vector{Matrix{Float64}}(undef, Nspin)
    UrotC = Vector{Matrix{Float64}}(undef, Nspin)
    for ikspin in 1:Nspin
        rotPrev[ikspin] = Matrix(1.0*I(Nstates))
        rotPrevC[ikspin] = Matrix(1.0*I(Nstates))
        rotPrevCinv[ikspin] = Matrix(1.0*I(Nstates))
        Urot[ikspin] = Matrix(1.0*I(Nstates))
        UrotC[ikspin] = Matrix(1.0*I(Nstates))
    end
    return RotationsCache(
        rotPrev,
        rotPrevC,
        rotPrevCinv,
        Urot,
        UrotC,
    )
end

function transform_psi_Haux_update_ebands!(
    Ham, psis, Haux, rots_cache
)
    Nspin = Ham.electrons.Nspin
    ebands = Ham.electrons.ebands
    hx = Ham.grid.hx
    #
    Urot = rots_cache.Urot
    UrotC = rots_cache.UrotC
    rotPrev = rots_cache.rotPrev
    rotPrevC = rots_cache.rotPrevC
    rotPrevCinv = rots_cache.rotPrevCinv
    #
    for ispin in 1:Nspin
        # ebands will be updated here
        ebands[:,ispin], Urot[ispin] = eigen(Hermitian(Haux[ispin]))
        # overwrite_Haux
        Haux[ispin] = diagm( 0 => Ham.electrons.ebands[:,ispin] )
        UrotC[ispin] = inv(sqrt(psis[ispin]' * psis[ispin])) ./ sqrt(hx)
        UrotC[ispin] = UrotC[ispin]*Urot[ispin] # extra rotation
        psis[ispin] = psis[ispin]*UrotC[ispin]
    end
    #
    for ispin in 1:Nspin
        # Save previous
        rotPrev[ispin] = rotPrev[ispin] * Urot[ispin]
        rotPrevC[ispin] = rotPrevC[ispin] * UrotC[ispin]
        rotPrevCinv[ispin] = inv(UrotC[ispin]) * rotPrevCinv[ispin]
    end
    return
end




function transform_psi_Haux!(psi, Haux)
    λ, Urot = eigen(Hermitian(Haux))
    psi[:,:] = psi[:,:]*Urot
    Haux[:,:] = diagm(0 => λ)
    return Urot
end

# Just like transform_psi_Haux + orthonormalization
function prepare_psi_Haux!(psi, Haux, hx)
    Udagger = inv(sqrt(psi'*psi)) ./ sqrt(hx)
    psi[:,:] = psi*Udagger
    # Make Haux diagonal
    λ, Urot = eigen(Hermitian(Haux))
    Haux[:,:] = diagm(0 => λ)
    psi[:,:] = psi[:,:]*Urot
    return Udagger, Urot # need these?
end

# This is reverse operation to prepare_psi_Haux
# XXX Probably not needed?
function revert_psi_Haux!(psi, Haux, Udagger, Urot)
    psi[:,:] = psi[:,:]*Urot'
    Haux[:,:] = Urot * Haux * Urot'
    psi[:,:] = psi * inv(Udagger)
    return
end



# Various functions to update Hamiltonian
#
# Input: ebands
# Modifies: ebands (copy), Focc, E_fermi, mTS
function update_from_ebands!(Ham, ebands)
    # Set ebands
    Ham.electrons.ebands[:] .= ebands[:]
    update_from_ebands!(Ham)
    return
end

# Update occupation numbers, calculate E_fermi and mTS
# ebands is read from Ham
function update_from_ebands!(Ham)
    # Set ebands
    ebands = Ham.electrons.ebands
    # Update occupation numbers
    kT = Ham.electrons.kT
    Focc = Ham.electrons.Focc
    Nelectrons = Ham.electrons.Nelectrons
    #
    Ham.electrons.E_fermi, Ham.energies.mTS = update_Focc!(
        Focc, smear_fermi, smear_fermi_entropy,
        ebands, Float64(Nelectrons), kT
    )
    return
end


# Input: psi
# Modifies: Rhoe, potentials, energies
function update_from_psis!(Ham, psis)

    Npoints = Ham.grid.Npoints
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    Vtot = Ham.potentials.Total
    rhoe = Ham.rhoe

    # Electron density
    calc_rhoe!(Ham, psis, rhoe)

    ispin = 1

    # Update the potentials
    # Note that Vxc, Vhartree, and Vtot refers to Ham.potentials
    ρ = reshape(rhoe, Npoints) # FIXME: need to do this is Poisson_solve_sum!
    # XXX probably use reduce
    Poisson_solve_sum!(Ham.grid, ρ, Vhartree)
    Vxc[:] = calc_Vxc_1d(Ham.xc_calc, rhoe)
    @views Vtot[:,ispin] = Vion[:] + Vhartree[:] + Vxc[:,ispin]

    return
end


# update_from_ebands! and update_from wavefunc!
# should be called before
# only Ham.energies is modified
# psi still need to be passed for calc_E_kin
function calc_Lfunc(
    Ham::Hamiltonian1d,
    psis::Vector{Matrix{Float64}}, # (Nbasis,Nstates)
)

    ebands = Ham.electrons.ebands
    # We don't support Nspin=2 yet
    @assert size(ebands,2) == 1

    ispin = 1
    psi = psis[ispin]

    mTS = Ham.energies.mTS
    hx = Ham.grid.hx    
    Vion = Ham.potentials.Ions
    Vhartree = Ham.potentials.Hartree
    rhoe = Ham.rhoe

    # Evaluate total energy components
    #
    # Kinetic energy
    Ekin = calc_E_kin(Ham, psi)
    Ham.energies.Kinetic = Ekin
    # Hartree
    Ehartree = 0.5*dot(rhoe[:,ispin], Vhartree)*hx
    Ham.energies.Hartree = Ehartree
    # Ion
    Eion = dot(rhoe, Vion)*hx
    Ham.energies.Ion = Eion
    # XC
    epsxc = calc_epsxc_1d(Ham.xc_calc, rhoe[:,ispin])
    Exc = dot(rhoe, epsxc)*hx
    Ham.energies.XC = Exc
    #
    # The total energy (also include nuclei-nuclei or ion-ion energy)
    Etot = Ekin + Ehartree + Eion + Exc + Ham.energies.NN + mTS

    return Etot
end


# The inputs are:
# - wavefunction psi, and
# - auxiliary Hamiltonian in diagonal form, stored as matrix with size (Nstates,Nspin)
#   Nspin=1 is assumed for the moment.
# Some fields of Ham will be modified
function calc_Lfunc_ebands!(
    Ham::Hamiltonian1d,
    psi::Matrix{Float64}, # (Nbasis,Nstates)
    ebands::Matrix{Float64} # (Nstates,Nspin)
)
    update_from_ebands!(Ham, ebands)
    update_from_wavefunc!(Ham, psi)
    Etot = calc_Lfunc(Ham, psi)
    return Etot
end


# The inputs are:
# - wavefunction psi, and
# - auxiliary Hamiltonian Haux. No support for spin polarization for now.
#
# Some fields of Ham will be modified
#
# psi and Haux must be transformed simultaneously by using some unitary matrix.
# The transformation chosen such that Haux transformed to diagonal form using
# eigendecomposition.
#
# psi and Haux will not be modified in place upon calling this function.
function calc_Lfunc_Haux!(
    Ham::Hamiltonian1d,
    psi, # (Nbasis,Nstates)
    Haux # (Nstates,Nstates)
)
    # Calculate ebands first
    ebands = zeros(Float64, size(psi,2), 1) # ebands need to be of size (Nstates,1)
    # Ham.electrons.ebands also can be used
    ebands[:,1], Urot = eigen(Hermitian(Haux)) # Force Haux to be Hermitian

    return calc_Lfunc_ebands!(Ham, psi*Urot, ebands)
end