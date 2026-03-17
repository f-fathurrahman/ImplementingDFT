# Taken from PWDFT.jl
# Some methods and references to PWDFT.Atoms and pspots are removed

"""
The type for describing electronic variables such as number of
electrons, occupation numbers, and energy levels.
"""
mutable struct Electrons
    Nelectrons::Float64
    Nstates::Int64
    Nstates_occ::Int64
    Focc::Array{Float64,2}
    ebands::Array{Float64,2}
    Nspin_wf::Int64
    use_smearing::Bool
    kT::Float64
    noncollinear::Bool
    E_fermi::Float64
    Nspin_dens::Int64
    domag::Bool
end

"""
Creates a 'dummy' instance of `Electrons` with only one electron.
"""
function Electrons()
    Nelectrons = 1
    Nstates = 1
    Nstates_occ = 1
    Focc = zeros(Nstates,1) # Nkpt=1
    ebands = zeros(Nstates,1) # use Nkpt=1
    Nspin_wf = 1
    Nspin_dens = 1
    use_smearing = false
    kT = 0.0
    noncollinear = false
    E_fermi = 0.0
    domag = false
    return Electrons(
        Nelectrons, Nstates, Nstates_occ, Focc, ebands, Nspin_wf,
        use_smearing, kT, noncollinear, E_fermi, Nspin_dens, domag
    )
end


function calc_Nstates(
    Nelectrons;
    Nstates_empty = nothing,
    use_smearing = false,
    noncollinear = false
)
    if noncollinear
        degspin = 1
    else
        degspin = 2
    end
    #
    Nstates = ceil(Int64, Nelectrons/degspin)
    if use_smearing && isnothing(Nstates_empty)
        # metallic case: add 20% more bands, with a minimum of 4
        Nstates = maximum([
            ceil(Int64, 1.2*Nelectrons/degspin),
            Nstates + 4
        ])
    else
        if Nstates_empty > 0
            Nstates += Nstates_empty
        end
    end
    #
    if isnothing(Nstates_empty)
        Nstates_empty = Nstates - ceil(Int64, Nelectrons/degspin)
        @assert Nstates_empty >= 0
    end
    @assert !isnothing(Nstates_empty)
    return Nstates, Nstates_empty
end


# Docs removed

function Electrons(
    atoms::Atoms1d;
    Nspin_wf = 1,
    Nkpt = 1,
    Nstates = nothing,
    Nstates_empty = nothing,
    noncollinear = false,
    domag = false,
    use_smearing = false
)
    if !noncollinear
        @assert Nspin_wf <= 2
    end

    # Determine Nspin_dens
    Nspin_dens = 1 # default
    # Collinear magnetism
    if !noncollinear && (Nspin_wf == 2)
        Nspin_dens = 2
    end
    # For noncollinear Nspin_wf = 1 and Nspin_dens = 4
    if noncollinear
        @assert Nspin_wf == 1
        Nspin_dens = 4
    end

    Nelectrons = get_Nelectrons(atoms)

    is_odd = round(Int64,Nelectrons)%2 == 1

    #XXX This is quite convoluted

    if isnothing(Nstates) && use_smearing
        # automatic determination of Nstates and Nstates_empty in case of smearing
        Nstates, Nstates_empty = calc_Nstates(
            Nelectrons,
            Nstates_empty = Nstates_empty,
            use_smearing = use_smearing,
            noncollinear = noncollinear
        )
    else
        # XXX The original logic, Nstates and Nstates_empty are
        if !noncollinear
            if isnothing(Nstates)
                Nstates = round(Int64, Nelectrons/2)
                if Nstates*2 < Nelectrons
                    Nstates = Nstates + 1
                end
                if !isnothing(Nstates_empty)
                    Nstates = Nstates + Nstates_empty
                else
                    Nstates_empty = 0 # Given Nstates_empty a valid value
                end
            else
                # Nstates is not given its default (invalid value)
                if !isnothing(Nstates_empty)
                    # In this case Nstates and Nstates_empty are both given
                    # We only allow to specify either Nstates or Nstates_empty
                    # We will error in this case
                    error("Please specify Nstates only or Nstates_empty only")
                end
                NstatesMin = round(Int64, Nelectrons/2)
                if NstatesMin*2 < Nelectrons
                    NstatesMin += 1
                end
                if Nstates < NstatesMin
                    error("Given Nstates is not enough: Nstates = $(Nstates), NstatesMin = $(NstatesMin)")
                end
                if Nstates > NstatesMin
                    Nstates_empty = Nstates - NstatesMin
                else
                    Nstates_empty = 0
                end
            end
            # Nstates, Nstates_empty must be set up to their valid values now
        else
            if isnothing(Nstates)
                error("Please specify Nstates manually")
            end
            # Nstates is assumed to be given
            Nstates_empty = Int64(Nstates - Nelectrons)
            @assert Nstates_empty >= 0 # cannot have negative Nstates_empty
        end
    end


    Focc = zeros(Float64, Nstates, Nkpt*Nspin_wf)
    ebands = zeros(Float64, Nstates, Nkpt*Nspin_wf)
    Nstates_occ = Nstates - Nstates_empty

    @info "Nstates = $(Nstates)"
    @info "Nstates_empty = $(Nstates_empty)"
    @info "Nstates_occ = $(Nstates_occ)"
    @info "Nspin_wf = $(Nspin_wf)"
    @info "Nspin_dens = $(Nspin_dens)"

    if noncollinear
        OCC_MAX = 1.0
    else
        if Nspin_wf == 1
            OCC_MAX = 2.0
        else
            OCC_MAX = 1.0
        end
    end
    @info "OCC_MAX = $(OCC_MAX)"

    if Nspin_wf == 1
        for ist in 1:Nstates_occ-1
            Focc[ist,:] .= OCC_MAX
        end
        if is_odd
            Focc[Nstates_occ,:] .= 1.0
        else
            Focc[Nstates_occ,:] .= OCC_MAX
        end
    else
        for ist = 1:Nstates_occ-1
            Focc[ist,:] .= OCC_MAX
        end
        idx_up = 1:Nkpt
        if is_odd
            # assign only to the spin up
            Focc[Nstates_occ,idx_up] .= OCC_MAX
        else
            Focc[Nstates_occ,:] .= OCC_MAX
        end
    end

    _check_Focc(Focc, Nkpt, Nelectrons)

    kT = 0.0
    E_fermi = 0.0
    return Electrons(
        Nelectrons, Nstates, Nstates_occ, Focc, ebands, Nspin_wf,
        use_smearing, kT, noncollinear, E_fermi, Nspin_dens,
        domag
    )
end


function _check_Focc(Focc::Matrix{Float64}, Nkpt::Int64, Nelectrons)
    sFocc = sum(Focc)/Nkpt
    # Check if the generated Focc is consistent
    if abs( sFocc - Nelectrons ) > eps()
        error("sum(Focc) = $sFocc and Nelectrons = $Nelectrons does not match")
    end
    return
end


function get_Nelectrons( atoms::Atoms1d )
    return sum(atoms.Zvals)
end
