mutable struct Electrons
    Nelectrons::Int64
    Nstates_occ::Int64
    Nstates_extra::Int64
    Nstates::Int64
    Nspin::Int64
    Focc::Matrix{Float64}
    ebands::Matrix{Float64}
    E_fermi::Float64
    kT::Float64
end


function Electrons(atoms::Atoms1d; Nstates_extra=0, Nspin=1)
    Nelectrons = sum(atoms.Zvals)
    Nstates_occ = round(Int64, Nelectrons/2)
    Nstates = Nstates_occ + Nstates_extra
    #
    is_odd = round(Int64,Nelectrons)%2 == 1
    
    Focc = zeros(Float64, Nstates, Nspin)
    if Nspin == 1
        # XXX Broadcast is used here in case we have second dimension of Nkpt*Nspin
        for ist in 1:Nstates_occ-1
            Focc[ist,:] .= 2.0
        end
        if is_odd
            Focc[Nstates_occ,:] .= 1.0
        else
            Focc[Nstates_occ,:] .= 2.0
        end
    end

    #
    ebands = zeros(Float64, Nstates, Nspin)

    E_fermi = 0.0
    kT = 0.0
    return Electrons(
        Nelectrons, Nstates_occ, Nstates_extra, Nstates,
        Nspin,
        Focc, ebands, E_fermi, kT
    )
end