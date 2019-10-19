function calc_rhoe( Ham::Hamiltonian, psi::Array{Float64,2} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)

    Rhoe = zeros(Float64,Nbasis)
    Focc = Ham.electrons.Focc

    for ist in 1:Nstates
        for ip in 1:Nbasis
            Rhoe[ip] = Rhoe[ip] + Focc[ist]*psi[ip,ist]*psi[ip,ist]
        end
    end
    return Rhoe
end
