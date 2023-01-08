mutable struct Energies
    Kinetic::Float64
    Ion::Float64
    Hartree::Float64
    XC::Float64
    NN::Float64
    mTS::Float64
end

function Energies()
    return Energies(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end


import Base: show
function show( io::IO, energies::Energies )

    @printf(io, "Kinetic    energy: %18.10f\n", energies.Kinetic )
    @printf(io, "Ion        energy: %18.10f\n", energies.Ion )
    @printf(io, "Hartree    energy: %18.10f\n", energies.Hartree )
    @printf(io, "XC         energy: %18.10f\n", energies.XC )
    @printf(io, "-TS              : %18.10f\n", energies.mTS)

    @printf(io, "-------------------------------------\n")
    
    E_elec = energies.Kinetic + energies.Ion +
             energies.Hartree + energies.XC + energies.mTS
    
    @printf(io, "Electronic energy: %18.10f\n", E_elec)
    @printf(io, "NN         energy: %18.10f\n", energies.NN )
    @printf(io, "-------------------------------------\n")
    
    E_total = E_elec + energies.NN
    
    @printf(io, "Total free energy: %18.10f\n", E_total)
    @printf(io, "\n")
    @printf(io, "Total energy (extrapolated to T=0): %18.10f\n", E_total - 0.5*energies.mTS)
    return
end

show( energies::Energies ) = show( stdout, energies )

import Base: println
println( energies::Energies) = show( energies )