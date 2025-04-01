function calc_rhoe( Ham::Hamiltonian1d, psi::Array{Float64,2} )
    Npoints = size(psi, 1)
    @assert Ham.electrons.Nspin == 1
    Nspin = 1 # XXX
    rhoe = zeros(Float64,Npoints,Nspin)
    calc_rhoe!(Ham, psi, rhoe)
    return rhoe
end

# In-place version
function calc_rhoe!(
    Ham::Hamiltonian1d,
    psi::Array{Float64,2},
    rhoe::Array{Float64,2}
)
    Npoints = size(psi, 1)
    Nstates = size(psi, 2)
    fill!(rhoe, 0.0)
    Focc = Ham.electrons.Focc
    @assert Ham.electrons.Nspin == 1
    ispin = 1 # XXX read ispin from Ham ?
    for ist in 1:Nstates, ip in 1:Npoints
        rhoe[ip,ispin] += Focc[ist,ispin] * psi[ip,ist] * psi[ip,ist]
    end
    return
end

# In-place version
function calc_rhoe!(
    Ham::Hamiltonian1d,
    psis::Vector{Matrix{Float64}},
    rhoe::Array{Float64,2}
)
    Nspin = size(psis, 1)
    @assert Nspin == Ham.electrons.Nspin
    Npoints = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    fill!(rhoe, 0.0)
    for ispin in 1:Nspin
        psi = psis[ispin]
        for ist in 1:Nstates, ip in 1:Npoints
            rhoe[ip,ispin] += Focc[ist,ispin] * psi[ip,ist] * psi[ip,ist]
        end
    end
    return
end