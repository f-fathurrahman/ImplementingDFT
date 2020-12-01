push!(LOAD_PATH, "./")

using ModelSCF1d

function main()
    # getting all the parameters
    dx = 0.5
    Nunit = 8   # number of units
    Lat = 10     # size of the lattice
    Ls = Nunit*Lat
    Ns = round(Integer, Ls / dx) # number of discretization points

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

    println(atoms)
end

main()