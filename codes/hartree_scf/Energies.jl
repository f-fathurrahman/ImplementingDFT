mutable struct Energies
    Kinetic::Float64
    Ps_loc::Float64
    Hartree::Float64
    NN::Float64
end

function Energies()
    return Energies(0.0, 0.0, 0.0, 0.0)
end

import Base: sum
function sum( ene::Energies )
    return ene.Kinetic + ene.Ps_loc + ene.Hartree + ene.NN
end

import Base: println
function println( ene::Energies )
    @printf("----------------------------\n")
    @printf("Total energy components\n")
    @printf("----------------------------\n")
    @printf("Kinetic = %18.10f\n", ene.Kinetic)
    @printf("Ps_loc  = %18.10f\n", ene.Ps_loc)
    @printf("Hartree = %18.10f\n", ene.Hartree)
    @printf("NN      = %18.10f\n", ene.NN)
    @printf("----------------------------\n")
    @printf("Total   = %18.10f\n", sum(ene))
end