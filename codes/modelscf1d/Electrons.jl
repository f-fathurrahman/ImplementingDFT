mutable struct Electrons
    Nelectrons::Int64
    Nstates_occ::Int64
    Nstates_extra::Int64
    Nstates::Int64
end


function Electrons(atoms::Atoms1d; Nstates_extra=0)
    Nelectrons = sum(atoms.Zvals)
    Nstates_occ = round(Int64, Nelectrons/2)
    Nstates = Nstates_occ + Nstates_extra
    return Electrons(Nelectrons, Nstates_occ, Nstates_extra, Nstates)
end