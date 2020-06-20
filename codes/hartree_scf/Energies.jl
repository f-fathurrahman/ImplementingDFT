mutable struct Energies
    Kinetic::Float64
    Ps_loc::Float64
    Hartree::Float64
end

function Energies()
    return Energies(0.0, 0.0, 0.0)
end

import Base: sum
function sum( ene::Energies )
    return ene.Kinetic + ene.Ps_loc + ene.Hartree
end