mutable struct Electrons
    Nelectrons::Int64
    Nstates::Int64
    Nstates_occ::Int64
    Focc::Array{Float64,1}
    ene::Array{Float64,1}
end

function Electrons( Nelectrons::Int64; Nstates_extra=0 )

    is_odd = (Nelectrons%2 == 1)

    Nstates_occ = round(Int64, Nelectrons/2)
    if is_odd
        Nstates_occ = Nstates_occ + 1
    end
    
    Nstates = Nstates_occ + Nstates_extra
    
    Focc = zeros(Float64,Nstates)
    ene = zeros(Float64,Nstates)

    if !is_odd
        for i in 1:Nstates_occ
           Focc[i] = 2.0 
        end
    else
        for i in 1:Nstates_occ-1
            Focc[1] = 2.0
        end
        Focc[Nstates_occ] = 1.0
    end

    return Electrons(Nelectrons, Nstates, Nstates_occ, Focc, ene)
end