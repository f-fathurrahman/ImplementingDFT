include("update_Hamiltonian.jl")

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

    # We don't support Nspin=2 yet
    @assert size(ebands,2) == 1

    update_from_ebands!(Ham, ebands)
    update_from_wavefunc!(Ham, psi)

    hx = Ham.grid.hx    
    Npoints = Ham.grid.Npoints
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    Vtot = Ham.potentials.Total
    rhoe = Ham.rhoe
    mTS = Ham.energies.mTS

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