function get_occ(ev::Vector{Float64}, nocc::Int64, Tbeta::Float64)
    #
    # usage: [occ,efermi] = get_occ(ev, nocc, Tbeta);
    #
    # Obtained from KSSOLV
    #
    # ZX: Transferred to julia language.
    nev = length(ev)
    tol = 1e-15
    maxiter = 100

    if (nev > nocc)
       #
       # use bisection to find efermi such that
       #       sum_i fermidirac(ev(i)) = nocc
       #
       ilb = nocc-1
       iub = nocc+1
       lb = ev[ilb]
       ub = ev[iub]
       #
       # make sure flb < nocc and fub > nocc
       #
       flb = sum(fermidirac(ev,lb,Tbeta))
       fub = sum(fermidirac(ev,ub,Tbeta))
       while ( (nocc-flb)*(fub-nocc) < 0 )
          #fprintf('getocc: initial bounds are off:\n');
          #fprintf('flb = #11.3e, fub = #11.3e, nocc = #d\n', flb,fub,nocc);
          if (flb > nocc)
             if (ilb > 1)
                ilb = ilb - 1
                lb = ev[ilb]
                flb = sum(fermidirac(ev,lb,Tbeta))
             else
                error("getocc: cannot find a lower bound for efermi, something is wrong")
                break
             end
          end;
          if (fub < nocc)
             if (iub < nev)
                iub = iub + 1
                ub  = ev[iub]
                fub = sum(fermidirac(ev,ub,Tbeta))
             else
                error("getocc: cannot find an upper bound for efermi, something is wrong, try increasing the number of wavefunctions in X0")
                break
             end
          end
       end
       #fprintf('flb = #11.3e, fub = #11.3e\n', flb, fub);
       #pause;
       efermi = (lb+ub)/2
       occ = fermidirac(ev,efermi,Tbeta)
       occsum = sum(occ)
       iter = 1
       while ( abs(occsum-nocc) > tol ) && ( iter < maxiter )
          if (occsum < nocc)
             lb = efermi
          else
             ub = efermi
          end
          efermi = (lb+ub)/2
          occ = fermidirac(ev,efermi,Tbeta)
          occsum = sum(occ)
          iter = iter + 1
       end
    elseif (nev == nocc)
       occ    = ones(nocc,1)
       efermi = ev[nocc]
    else
       occ = []
       efermi = 0
       error("The number of eigenvalues in ev should be larger than nocc")
    end
    return occ, efermi
end

function fermidirac(ev, efermi, Tbeta::Float64)
    f = 1 ./(1 .+exp.(Tbeta*(ev .- efermi)))
    return f
end
