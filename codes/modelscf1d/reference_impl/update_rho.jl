function update_rho!(H::Hamiltonian, nocc::Int64)

    ev = H.ev
    (occ, fermi) = get_occ(ev, nocc, H.Tbeta)

    occ = occ * H.nspin
    rho = sum( H.psi.^2 * Matrix(Diagonal(occ)), dims=2 ) / H.dx
    println("integ rho = ", sum(rho)*H.dx)

    # Total energy
    E = sum(ev .* occ)

    # Helmholtz free energy
    intg = H.Tbeta*(fermi .- ev)
    ff = zeros(H.Neigs,1)
    for i = 1 : H.Neigs
      if( intg[i] > 30 )  # Avoid numerical problem.
        ff[i] = ev[i]-fermi
      else
        ff[i] = -1/H.Tbeta * log(1.0 + exp(intg[i]))
      end
    end
    F = sum(ff.*occ) + fermi * nocc * H.nspin


    println("\nEigenvalues and occupations")
    println("nocc = ", nocc)
    for ist in 1:H.Neigs
        @printf("%3d %18.10f %18.10f\n", ist, ev[ist], occ[ist])
    end
    println("Fermi energy = ", fermi)
    println("Free energy  = ", F)
    exit()


    H.occ = occ
    H.fermi = fermi
    H.Eband = E
    H.Fband = F
    H.rho = rho
end

function update_rhoa!(H::Hamiltonian)
    H.rhoa, H.drhoa = pseudocharge(H.gridpos, H.Ls, H.atoms, H.YukawaK, H.epsil0)
end