mutable struct Electrons
    Nelectrons::Int64
    Nstates::Int64
    Nstates_occ::Int64
    Focc::Array{Float64,1}
    eorbs::Array{Float64,1}
end


function Electrons(
    atoms::Atoms, pspots::Array{PsPot_GTH,1};
    Nstates=nothing, Nstates_empty=0
)

    Nelectrons = get_Nelectrons(atoms,pspots)

    is_odd = round(Int64,Nelectrons)%2 == 1

    # If Nstates is not specified and Nstates_empty == 0, we calculate
    # Nstates manually from Nelectrons
    if (Nstates == nothing)
        Nstates = round( Int64, Nelectrons/2 )
        if Nstates*2 < Nelectrons
            Nstates = Nstates + 1
        end
        if Nstates_empty > 0
            Nstates = Nstates + Nstates_empty
        end
    end

    Focc = zeros(Float64,Nstates)
    eorbs = zeros(Float64,Nstates)
    
    Nstates_occ = Nstates - Nstates_empty
    
    for ist = 1:Nstates_occ-1
        Focc[ist] = 2.0
    end
    if is_odd
        Focc[Nstates_occ] = 1.0
    else
        Focc[Nstates_occ] = 2.0
    end

    sFocc = sum(Focc)
    # Check if the generated Focc is consistent
    if abs( sFocc - Nelectrons ) > eps()
        @printf("sum Focc = %f, Nelectrons = %f\n", sFocc, Nelectrons)
        error(@sprintf("ERROR: diff sum(Focc) and Nelectrons is not small\n"))
    end

    return Electrons( Nelectrons, Nstates, Nstates_occ, Focc, eorbs )
end
