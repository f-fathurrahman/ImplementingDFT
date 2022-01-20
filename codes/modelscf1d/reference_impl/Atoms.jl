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
     function Atoms(Natoms, R, sigma, omega, Eqdist, mass, Z, nocc, force)
         return new(Natoms, R, sigma, omega, Eqdist, mass, Z, nocc, force)
     end
     #TODO:   this is a hack to make this thing work we need to come back and find
     #         how to properly define this function
     function Atoms(Natoms, R, sigma, omega, Eqdist, mass, Z, nocc)
         return new(Natoms, R, sigma, omega, Eqdist, mass, Z, nocc, 0*nocc)
     end
end
