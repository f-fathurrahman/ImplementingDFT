function calc_rhoe!(
    Ham::Hamiltonian1d,
    psi::Array{Float64,2},
    rhoe::Array{Float64,2}
)
    Npoints = size(psi, 1)
    Nstates = size(psi, 2)
    fill!(rhoe, 0.0)
    Focc = Ham.electrons.Focc
    for ist in 1:Nstates, ip in 1:Npoints
        rhoe[ip] += Focc[ist] * psi[ip,ist] * psi[ip,ist]
    end
    return
end
