function calc_rhoe( Ham::Hamiltonian, psi::Array{Float64,2} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    Focc = Ham.electrons.Focc
    Rhoe = zeros(Float64,Nbasis)
    for ist in 1:Nstates
        for ip in 1:Nbasis
            Rhoe[ip] = Rhoe[ip] + Focc[ist]*psi[ip,ist]*psi[ip,ist]
        end
    end
    return Rhoe
end


function calc_rhoe!( Ham::Hamiltonian, psi::Array{Float64,2}, Rhoe::Array{Float64,1} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    Focc = Ham.electrons.Focc
    Rhoe[:] = zeros(Float64,Nbasis)
    for ist in 1:Nstates
        for ip in 1:Nbasis
            Rhoe[ip] = Rhoe[ip] + Focc[ist]*psi[ip,ist]*psi[ip,ist]
        end
    end
    return
end