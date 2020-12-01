mutable struct Atoms
     # atoms struct to store the data from the atoms
     Natoms::Int64
     R::Array{Float64,2}
     sigma::Array{Float64,2}
     omega::Array{Float64,2}
     Eqdist::Array{Float64,2}
     mass::Array{Float64,2}
     Z::Array{Float64,2}
     nocc::Array{Int64,2}
     force::Array{Float64,2}
end

function Atoms(Natoms, R, sigma, omega, Eqdist, mass, Z, nocc)
    return Atoms(Natoms, R, sigma, omega, Eqdist, mass, Z, nocc, 0*nocc)
end