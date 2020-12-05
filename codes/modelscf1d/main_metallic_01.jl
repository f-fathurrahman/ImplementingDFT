push!(LOAD_PATH, "./")

using Printf
using ModelSCF1d

function main()
    dx = 0.5
    Nunit = 8   # number of units
    Lat = 10     # size of the lattice
    Ls = Nunit*Lat
    # using the default values in Lin's code
    YukawaK = 0.0100
    n_extra = 50
    epsil0 = 10.0
    T_elec = 50000.0

    kb = 3.1668e-6
    au2K = 315774.67
    Tbeta = au2K / T_elec

    betamix = 0.4
    mixdim = 10

    Ndist  = 1   # Temporary variable
    Natoms = round(Integer, Nunit / Ndist) # number of atoms

    sigma  = ones(Natoms,1)*(3.0)  # metal
    omega  = ones(Natoms,1)*0.03
    Eqdist = ones(Natoms,1)*10.0
    mass   = ones(Natoms,1)*42000.0
    nocc   = ones(Natoms,1)*2          # number of electrons per atom
    Z      = nocc

    x0 = zeros(Natoms) # this is defined as an 2D array
    for j in 1:Natoms
        x0[j] = (j - 0.5)*Lat*Ndist + dx
    end

    R = reshape(x0, length(x0), 1) # we need to make it a two dimensional array
    # creating an atom structure
    atoms = Atoms(Natoms, R, sigma,  omega,  Eqdist, mass, Z, nocc)
    # allocating a Hamiltonian
    Ham = Hamiltonian(Lat, Nunit, n_extra, dx, atoms, YukawaK, epsil0, Tbeta)

    # total number of occupied orbitals
    Nocc = round(Integer, sum(atoms.nocc) / Ham.nspin)

    # setting the options for the scf iteration
    KerkerB = 0.5
    mixOpts = AndersonMixOptions(Ham.Ns, betamix, mixdim)
    eigOpts = EigensolverOptions(1.e-8, 1000, "eigs")
    scfOpts = SCFOptions(1.e-7, 100, eigOpts, mixOpts)

    # initialize the potentials within the Hemiltonian, setting H[\rho_0]
    init_pot!(Ham, Nocc)

    # running the scf iteration
    VtoterrHist = scf_potmix!(Ham, scfOpts)
    
    println("bandgap = ", Ham.ev[Nocc+1] - Ham.ev[Nocc])
end

main()