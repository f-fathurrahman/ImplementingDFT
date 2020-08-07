function update_Focc!(
    Focc::Array{Float64,1},
    smear_func, smear_func_entropy,
    evals::Array{Float64,1},
    Nelectrons::Float64,
    kT::Float64,
)

    E_f = find_E_fermi( smear_func, evals, Nelectrons, kT )

    Nstates = size(evals,1)
    w = 2.0 # weight factor
    for ist in 1:Nstates
        Focc[ist] = w*smear_func( evals[ist], E_f, kT )
    end
  
    mTS = 0.0
    for ist = 1:Nstates
        mTS = mTS - w*kT*smear_func_entropy( evals[ist], E_f, kT )
    end
    return E_f, mTS  
end


function sum_Focc(
    smear_func,
    evals::Array{Float64,1},
    ene::Float64,
    kT::Float64,
)
    Nstates = size(evals,1)
    ss = 0.0
    for ist = 1:Nstates
        ss = ss + smear_func( evals[ist], ene, kT )
    end
    return 2.0*ss # weight=2.0 (Nspin=1, assuming doubly occupied)
end


function find_E_fermi(
    smear_func,
    evals::Array{Float64,1},
    Nelectrons::Float64,
    kT::Float64;
    NiterMax=300, verbose=false
)

    Nstates = size(evals,1)

    # determine lower and upper bound for bisection
    Elw = evals[1]
    Eup = evals[Nstates]
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