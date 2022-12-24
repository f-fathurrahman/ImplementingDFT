mutable struct Electrons
    Nelectrons::Int64
    Nstates_occ::Int64
    Nstates_extra::Int64
    Nstates::Int64
    Focc::Matrix{Float64}
    ebands::Matrix{Float64}
end


function Electrons(atoms::Atoms1d; Nstates_extra=0, Nspin=1)
    Nelectrons = sum(atoms.Zvals)
    Nstates_occ = round(Int64, Nelectrons/2)
    Nstates = Nstates_occ + Nstates_extra
    #
    Focc = zeros(Float64, Nstates, Nspin)
    Focc[1:Nstates_occ] .= 2.0
    #
    ebands = zeros(Float64, Nstates, Nspin)
    #
    return Electrons(Nelectrons, Nstates_occ, Nstates_extra, Nstates, Focc, ebands)
end