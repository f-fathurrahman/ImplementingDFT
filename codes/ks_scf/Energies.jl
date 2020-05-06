mutable struct Energies
    Kinetic::Float64
    Ps_loc::Float64
    Hartree::Float64
    XC::Float64
end

function Energies()
    return Energies(0.0, 0.0, 0.0, 0.0)
end

import Base: sum
function sum( ene::Energies )
    return ene.Kinetic + ene.Ps_loc + ene.Hartree + ene.XC
end

import Base: println
function println( ene::Energies )

    @printf("----------------------------\n")
    @printf("Total energy components\n")
    @printf("----------------------------\n")
    @printf("Kinetic = %18.10f\n", ene.Kinetic)
    @printf("Ps_loc  = %18.10f\n", ene.Ps_loc)
    @printf("Hartree = %18.10f\n", ene.Hartree)
    @printf("XC      = %18.10f\n", ene.XC)
    @printf("----------------------------\n")
    @printf("Total   = %18.10f\n", sum(ene))
end