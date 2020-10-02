mutable struct Electrons
    Nspin::Int64
    Nelectrons::Int64
    Nstates::Int64
    Nstates_occ::Int64
    Focc::Array{Float64,2}
    eorbs::Array{Float64,2}
end


function Electrons(
    atoms::Atoms, pspots::Array{PsPot_GTH,1};
    Nstates=nothing, Nstates_extra=0, Nspin=1
)

    @assert Nspin <= 2

    Nelectrons = get_Nelectrons(atoms,pspots)

    is_odd = round(Int64,Nelectrons)%2 == 1

    # If Nstates is not specified and Nstates_extra == 0, we calculate
    # Nstates manually from Nelectrons
    if (Nstates == nothing)
        Nstates = round( Int64, Nelectrons/2 )
        if Nstates*2 < Nelectrons
            Nstates = Nstates + 1
        end
        if Nstates_extra > 0
            Nstates = Nstates + Nstates_extra
        end
    end

    Focc = zeros(Float64,Nstates,Nspin)
    eorbs = zeros(Float64,Nstates,Nspin)
    
    Nstates_occ = Nstates - Nstates_extra
    
    if Nspin == 1
        for ist = 1:Nstates_occ-1
            Focc[ist,1] = 2.0
        end
        if is_odd
            Focc[Nstates_occ,1] = 1.0
        else
            Focc[Nstates_occ,1] = 2.0
        end
    elseif Nspin == 2
        for ist = 1:Nstates_occ-1
            Focc[ist,1] = 1.0
            Focc[ist,2] = 1.0
        end
        if is_odd
            Focc[Nstates_occ,1] = 1.0
        end
    end

    sFocc = sum(Focc)
    # Check if the generated Focc is consistent
    if abs( sFocc - Nelectrons ) > eps()
        @printf("sum Focc = %f, Nelectrons = %f\n", sFocc, Nelectrons)
        error(@sprintf("ERROR: diff sum(Focc) and Nelectrons is not small\n"))
    end

    return Electrons( Nspin, Nelectrons, Nstates, Nstates_occ, Focc, eorbs )
end

import Base: show
function show( io::IO, electrons::Electrons; header=true, all_states=false )

    Nspin = electrons.Nspin
    Focc = electrons.Focc
    Nelectrons = electrons.Nelectrons
    Nstates = electrons.Nstates
    Nstates_occ = electrons.Nstates_occ
    Nspin = size(Focc,2)

    if header
        @printf(io, "\n")
        @printf(io, "                                    ---------\n")
        @printf(io, "                                    Electrons\n")
        @printf(io, "                                    ---------\n")
        @printf(io, "\n")
    end
    @printf(io, "Nspin         = %8d\n", Nspin)
    @printf(io, "Nelectrons    = %8d\n", Nelectrons)
    @printf(io, "Nstates       = %8d\n", Nstates)
    @printf(io, "Nstates_occ   = %8d\n", Nstates_occ)
    @printf(io, "Nstates_empty = %8d\n\n", Nstates - Nstates_occ)

    if Nspin == 1
        @printf(io, "Occupation numbers: (spin-paired)\n\n")
    else
        @printf(io, "Occupation numbers: (spin-polarized)\n\n")
    end

    if Nstates < 8
        all_states = true
    end

    Nk_per_line = 8

    if all_states
        for ist = 1:Nstates
            @printf(io, "state #%4d = ", ist)
            for ispin in 1:Nspin
                @printf(io, "%8.5f ", Focc[ist,ispin])
                if (ispin % Nk_per_line) == 0
                    @printf(io, "\n")
                    @printf(io, "              ")
                end
            end
            #
            @printf(io, "\n")
        end
    else
        for ist = 1:4
            @printf(io, "state #%4d = ", ist)
            for ispin in 1:Nspin
                @printf(io, "%8.5f ", Focc[ist,ispin])
                if (ispin % Nk_per_line) == 0
                    @printf(io, "\n")
                    @printf(io, "              ")
                end
            end
            @printf(io, "\n")
        end
        @printf(io, ".....\n\n")
        #
        for ist in Nstates-3:Nstates
            @printf(io, "state #%4d = ", ist)
            for ispin in 1:Nkspin
                @printf(io, "%8.5f ", Focc[ist,ispin])
                if (ispin % Nk_per_line) == 0
                    @printf(io, "\n")
                    @printf(io, "              ")
                end
            end
            @printf(io, "\n")
        end
    end
end

function show( electrons::Electrons; header=true, all_states=false )
    show( stdout, electrons, header=header, all_states=all_states )
end
