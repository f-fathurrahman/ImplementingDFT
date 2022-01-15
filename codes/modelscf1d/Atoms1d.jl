# atoms struct to store the data from the atoms
mutable struct Atoms1d
    Natoms::Int64
    positions::Vector{Float64}
    Zvals::Array{Float64}
    σ::Vector{Float64}
    masses::Vector{Float64}
    forces::Vector{Float64}
    L::Float64  # cell length
end

function Atoms1d( positions, Zvals, σ, masses, L )
    Natoms = length(positions)
    @assert Natoms == length(Zvals)
    @assert Natoms == length(σ)
    @assert Natoms == length(masses)
    forces = zeros(Natoms)
    #
    return Atoms1d( Natoms, positions, Zvals, σ, masses, forces, L )
end