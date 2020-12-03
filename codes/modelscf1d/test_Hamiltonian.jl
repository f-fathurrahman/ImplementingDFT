push!(LOAD_PATH, "./")

using Printf
using ModelSCF1d

function main()
    # getting all the parameters
    dx = 0.5
    Nunit = 8   # number of units
    Lat = 10     # size of the lattice
    Ls = Nunit*Lat
    Ns = round(Integer, Ls/dx) # number of discretization points

    Ndist = 1
    Natoms = round(Integer, Nunit / Ndist)
    
    sigma  = ones(Natoms,1)*(1.0)  # insulator
    omega  = ones(Natoms,1)*0.03
    Eqdist = ones(Natoms,1)*10.0
    mass   = ones(Natoms,1)*42000.0
    nocc   = ones(Natoms,1)*2     # number of electrons per atom
    Z      = nocc
    
    # defining the position of the nuclei
    R = reshape([ (j-0.5)*Lat*Ndist+dx for j=1:Natoms ], Natoms, 1)

    # creating an atom structure
    atoms = Atoms(Natoms, R, sigma, omega, Eqdist, mass, Z, nocc)

    # using the default values in Lin's code
    YukawaK = 0.0100
    n_extra = 10
    epsil0 = 10.0
    T_elec = 100.0

    kb = 3.1668e-6
    au2K = 315774.67
    Tbeta = au2K / T_elec

    betamix = 0.5
    mixdim = 10

    Ham = Hamiltonian(Lat, Nunit, n_extra, dx, atoms,YukawaK, epsil0, Tbeta)
    # total number of occupied orbitals
    Nocc = round(Integer, sum(atoms.nocc) / Ham.nspin)
    
    mixOpts = AndersonMixOptions(Ham.Ns, betamix, mixdim)
    eigOpts = EigensolverOptions(1.e-10, 1000, "eigs")
    #eigOpts = eigOptions(1.e-10, 1000, "lobpcg_sep")
    scfOpts = SCFOptions(1.e-8, 300, eigOpts, mixOpts)
    
    # initialize the potentials within the Hemiltonian, setting H[\rho_0]
    init_pot!(Ham, Nocc)

    println(typeof(Ham.occ))

    VtoterrHist = scf_potmix!(Ham, scfOpts)
    
    if VtoterrHist[end] > scfOpts.SCFtol
        println("convergence not achieved!! ")
    end
    
    println(length(VtoterrHist))

end

main()