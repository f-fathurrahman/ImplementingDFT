function calc_rhoe(
    Ham::Hamiltonian,
    psis::Vector{Array{Float64,2}}
)
    Nspin = size(psis, 1)
    Nbasis = size(psis[1], 1)
    Rhoe = zeros(Float64, Nbasis, Nspin)
    calc_rhoe!(Ham, psi, Rhoe)
    return Rhoe
end

function calc_rhoe!(
    Ham::Hamiltonian,
    psis::Vector{Array{Float64,2}},
    Rhoe::Array{Float64,2}
)
    Nspin = size(psis, 1)
    Nbasis = size(psis[1], 1)
    Nstates = size(psis[1], 2)
    Focc = Ham.electrons.Focc
    fill!(Rhoe, 0.0)
    for ispin in 1:Nspin
        psi = psis[ispin]
        for ist in 1:Nstates, ip in 1:Nbasis
            Rhoe[ip,ispin] = Rhoe[ip,ispin] + Focc[ist,ispin]*psi[ip,ist]*psi[ip,ist]
        end
    end
    return
end