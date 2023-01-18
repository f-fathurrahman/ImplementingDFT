struct OccupationUpdater
    smear_func::Function
    smear_func_entropy::Function
    smear_func_prime::Function
    kT::Float64
end

function OccupationUpdater(; smear_type=:FermiDirac, kT::Float64=0.01)
    if smear_type == :FermiDirac
        return OccupationUpdater(
            smear_fermi, smear_fermi_entropy, smear_fermi_prime, kT
        )
    elseif smear_type == :Gaussian
        return OccupationUpdater(
            smear_gauss, smear_gauss_entropy, smear_gauss_prime, kT
        )
    elseif (smear_type == :MarzariVanderbilt) || (smear_type == :cold)
        return OccupationUpdater(
            smear_cold, smear_cold_entropy, smear_cold_prime, kT
        )
    else
        println("Unknown smear_type: ", smear_type)
        error()
    end
end


function update_Focc!(
    occ_updater::OccupationUpdater,
    Focc::Matrix{Float64},
    evals::Matrix{Float64},
    Nelectrons
)
    return update_Focc!(Focc, occ_updater.smear_func, occ_updater.smear_func_entropy,
        Nelectrons, occ_updater.kT)
end

# Update occupation number according to `smear_func`
# with given energy eigenvalues `evals`
function update_Focc!(
    Focc::Array{Float64,2},
    smear_func, smear_func_entropy,
    evals::Array{Float64,2},
    Nelectrons,
    kT::Float64
)

    E_f = find_E_fermi( smear_func, evals, Nelectrons, kT )

    Nstates = size(evals,1)
    Nspin = size(evals,2)
    
    if Nspin == 1
        w = 2.0 # weight factor
    else
        w = 1.0
    end

    for ispin in 1:Nspin
        for ist in 1:Nstates
            Focc[ist,ispin] = w*smear_func( evals[ist,ispin], E_f, kT )
        end
    end
  
    mTS = 0.0
    for ispin in 1:Nspin
        for ist in 1:Nstates
            mTS = mTS - w*kT*smear_func_entropy( evals[ist,ispin], E_f, kT )
        end
    end
    return E_f, mTS  
end

# FIXME: include in update_Focc! ?
function calc_electronic_entropy( smear_func_entropy, evals, E_f, kT )
    Nstates = size(evals,1)
    Nspin = size(evals,2)
    mTS = 0.0
    if Nspin == 1
        w = 2.0 # weight factor
    else
        w = 1.0
    end
    for ispin in 1:Nspin
        for ist in 1:Nstates
            mTS = mTS - w*kT*smear_func_entropy( evals[ist,ispin], E_f, kT )
        end
    end
    return mTS
end


function sum_Focc(
    smear_func,
    evals::Array{Float64,2},
    ene::Float64,
    kT::Float64
)
    Nstates = size(evals,1)
    Nspin = size(evals,2)
    ss = 0.0
    for ispin in 1:Nspin
        for ist in 1:Nstates
            ss = ss + smear_func( evals[ist,ispin], ene, kT )
        end
    end
    if Nspin == 1
        return 2.0*ss
    else
        return ss
    end
end


function find_E_fermi(
    smear_func,
    evals::Array{Float64,2},
    Nelectrons,
    kT::Float64;
    NiterMax=300, verbose=false
)

    Nstates = size(evals,1)
    Nspin = size(evals,2)

    # determine lower and upper bound for bisection
    Elw = minimum(evals[1,:]) # minimum for all spin
    Eup = maximum(evals[Nstates,:]) # maximum for all spin

    Elw = Elw - 2*kT
    Eup = Eup + 2*kT

    verbose && println("Elw = ", Elw)
    verbose && println("Eup = ", Eup)

    sumlw = sum_Focc( smear_func, evals, Elw, kT )
    sumup = sum_Focc( smear_func, evals, Eup, kT )

    SMALL = 1e-10

    if ( (sumup - Nelectrons) < -eps() ) ||
       ( (sumlw - Nelectrons) >  eps() )
        @printf("sumup = %18.10f\n", sumup)
        @printf("sumlw = %18.10f\n", sumlw)
        @printf("Nelectrons = %18.10f\n", Nelectrons)
        error("Bounds for E_fermi is not found")
    end

    #
    # Start bisection
    #
    Ef = 0.5*(Eup + Elw)
    Ef_old = Ef
    for iter = 1:NiterMax
        sum_mid = sum_Focc( smear_func, evals, Ef, kT )
        if abs(sum_mid-Nelectrons) < SMALL
            return Ef
        elseif sum_mid-Nelectrons < -SMALL
            Elw = Ef
        else
            Eup = Ef
        end
        Ef = 0.5*(Eup + Elw)
        diff_Ef = abs(Ef-Ef_old)
        if verbose
            @printf("find_E_fermi: %3d %18.10f %18.10f %18.10e\n", iter, Ef, sum_mid, diff_Ef)
        end
        if diff_Ef < SMALL
            return Ef
        end
        Ef_old = Ef
    end

    @printf("WARNING: Ef is not found after %d iterations\n", NiterMax)
    return Ef
    
end