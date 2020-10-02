function update_Focc!(
    Focc::Array{Float64,2},
    smear_func, smear_func_entropy,
    evals::Array{Float64,2},
    Nelectrons::Float64,
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
    Nelectrons::Float64,
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