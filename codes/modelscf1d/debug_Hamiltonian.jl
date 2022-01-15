function debug_Hamiltonian(Lat, Nunit, n_extra, dx, atoms, YukawaK, epsil0, Tbeta)
    Ls = Nunit*Lat
    Ns = round(Int64, Ls/dx)
    dx = Ls/Ns

    println("Ls = ", Ls)
    println("Ns = ", Ns)
    println("dx = ", dx)

    # defining the grid
    gridpos = zeros(Ns) # allocating as a 2D Array
    for ip in 1:Ns
        gridpos[ip] = (ip-1)*dx
    end
    println(gridpos)

    posstart = 0
    posidx   = 0
    # Initialize the atom positions
    Neigs = sum(atoms.nocc) + n_extra
    println("Natoms = ", atoms.Natoms)
    println("n_extra = ", n_extra)
    println("Neigs = ", Neigs)
    println("nocc = ", atoms.nocc)

    Ns_glb = Ns

    # we define the Fourier multipliers as an 2D array
    Gx = zeros(Ns)
    # Follow the original ordering
    # TODO: Check the rounding
    ig = 0
    for i in 0:round(Int64,Ns/2-1)
        ig = ig + 1
        Gx[ig] = 2π/Ls * i
    end
    for i in round(Int64,-Ns/2):-1
        ig = ig + 1
        Gx[ig] = 2π/Ls * i
    end
    Gkin = 0.5 * Gx.^2
    
    println("ig = ", ig)
    println("Ns = ", Ns)

    for ig in 1:Ns
        @printf("%18.10f %18.10f\n", Gx[ig], Gkin[ig])
    end


    rhoa, drhoa = calc_pseudocharge(gridpos, Ls, atoms, YukawaK, epsil0)


    nspin = 1
    # From init_pot
    # nocc number of occupied states
    #function to initialize the potential in the Hamiltonian class
    rho  = -rhoa
    println("Integ rho = ", sum(rho)*dx)
    rho  = rho/(sum(rho)*dx) * (sum(atoms.nocc)*nspin) # renormalize
    println("Integ rho after = ", sum(rho)*dx)
    #Vhar = hartree_pot_bc(rho + rhoa, H)
    #Vtot = H.Vhar   # No exchange-correlation
    #Vtot = H.Vtot .- mean(H.Vtot)# IMPORTANT (zero mean?)


#=
    # TODO: we need to figure out the type of each of the fields to properlu
    # initialize them
    rho = zeros(1,1)
    Vhar = zeros(1,1)
    Vtot = zeros(1,1)
    # drhoa = []
    ev = []
    psi = zeros(1,1)
    fermi = 0.0
    occ = []
    Eband = 0.0
    Fband = 0.0
    Ftot = 0.0
=#


    return
end


function get_matrix( H::Hamiltonian )
    # create the matrix version of the Hmailtonian
    A = real( ifft(diagm(0 => H.kmul[:]) * fft(
            Matrix{Float64}(I, length(H.kmul), length(H.kmul)),1),1
        )
    )
    A += diagm(0 => H.Vtot[:])
    # we symmetrize A
    return 0.5*(A + A')
end