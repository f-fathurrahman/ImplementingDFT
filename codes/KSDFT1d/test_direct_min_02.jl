push!(LOAD_PATH, "./")

using Printf
using LinearAlgebra
using KSDFT1d

include("BroydenMixer.jl")

function create_atoms()
    Natoms = 3
    σ = ones(Float64, Natoms)*(1.0)
    masses = ones(Float64, Natoms)*42000.0
    Zvals = ones(Float64, Natoms)*5.0
    L = 10.0
    atpos = zeros(Float64, Natoms)
    atpos[1] = -1.5
    atpos[2] =  0.1
    atpos[3] =  1.5
    return Atoms1d( atpos, Zvals, σ, masses, L )
end

function pot_gaussian( x, x0 )
    return -25.0*exp(-4.5*(x-x0)^2)
end

function init_Vions!(Ham)
    atpos = Ham.atoms.positions
    Natoms = Ham.atoms.Natoms
    for ia in 1:Natoms
        Ham.potentials.Ions[:] += pot_gaussian.(Ham.grid.x, atpos[ia])
    end
    return
end

function init_Hamiltonian()
    atoms = create_atoms()
    Ham = Hamiltonian1d(atoms, 51, Nstates_extra=6)
    init_Vions!(Ham)
    return Ham
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
function calc_KohnSham_Etotal!(
    Ham,
    psi::Matrix{Float64},
    η::Matrix{Float64}
)

    kT = 0.1*eV2Ha # XXX FIXED

    # Calculate Kohn-Sham eigenvalues and occupation numbers
    Focc = Ham.electrons.Focc
    evals = Ham.electrons.ebands
    
    evals[:,1], rot = eigen(Hermitian(η)) # need to symmetrize?
    Nelectrons = Ham.electrons.Nelectrons
    E_f, mTS = update_Focc!( Focc, smear_fermi, smear_fermi_entropy,
                  evals, Float64(Nelectrons), kT )

    Ham.energies.mTS = mTS

    println("E_f = ", E_f)
    println("mTS = ", mTS)

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
    Etot = Ekin + Ehartree + Eion + Exc + Ham.energies.NN + mTS

    println(Ham.energies)

    return Etot
end

function prec_invK(Ham::Hamiltonian1d, v)
    return inv(Ham.Kmat)*v
end


#=
function main()

    Ham = init_Hamiltonian()

    hx = Ham.grid.hx
    Npoints = Ham.grid.Npoints
    Nelectrons = Ham.electrons.Nelectrons
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc

    Kmat = Ham.Kmat
    Vtot = Ham.potentials.Total
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    rhoe = Ham.rhoe

    Nspin = 1
    Hmat = zeros(Float64, Npoints, Npoints)
    epsxc = zeros(Float64, Npoints)
    rhoe_new = zeros(Float64, Npoints, Nspin)

    Etot = Inf
    Etot_old = Etot
    E_NN = calc_E_NN(Ham.atoms)

    Ekin = 0.0
    Ehartree = 0.0
    Eion = 0.0
    Exc = 0.0

    betamix = 0.5
    mixer = BroydenMixer(rhoe, betamix, mixdim=8)

    Focc = Ham.electrons.Focc
    use_smearing = true
    kT = 0.1*eV2Ha #0.1 # 0.3*eV2Ha
    mTS = 0.0

    println("kT = ", kT)

    evals = Ham.electrons.ebands


end

main()
=#