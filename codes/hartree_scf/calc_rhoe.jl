function calc_rhoe( Ham, psi::Array{Float64,2} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)

    Rhoe = zeros(Float64,Nbasis)

    for ist in 1:Nstates
        f = Ham.electrons.Focc[ist]
        for ip in 1:Nbasis
            Rhoe[ip] = Rhoe[ip] + f*psi[ip,ist]*psi[ip,ist]
        end
    end
    return Rhoe
end