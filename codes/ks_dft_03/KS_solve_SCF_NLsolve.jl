using NLsolve

"""
NLsolve driver.
Adapted from DFTK.jl
"""
function scf_NLsolve_solver( m=5; kwargs... )

    function fix_point_solver( f, x0, tol, NiterMax )
        res = NLsolve.nlsolve(x -> f(x) - x, x0;
            m=m, method=:anderson, xtol=tol, ftol=0.0, show_trace=true,
            iterations=NiterMax, kwargs... )
        ( fixpoint=res.zero, converged=NLsolve.converged(res) )
    end

    return fix_point_solver

end

function KS_solve_SCF_NLsolve!(
    Ham::Hamiltonian, psi::Array{Float64,2};
    NiterMax=200, betamix=0.5,
    etot_conv_thr=1e-6,
    diag_func=diag_LOBPCG!,
    kT=0.01
)

    Npoints = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol

    Rhoe_out = zeros(Float64,Npoints)
    Rhoe = zeros(Float64,Npoints)
    
    calc_rhoe!( Ham, psi, Rhoe )
    update!( Ham, Rhoe )
    
    evals = zeros(Float64,Nstates)
    Etot_old = 0.0
    dEtot = 0.0
    dRhoe = 0.0
    NiterMax = 100

    ethr_last=1e-5
    ethr = 0.1

    Ham.energies.NN = calc_E_NN( Ham.atoms, Ham.pspots )
    Etot_old = 0.0

    function density_map( Rhoe_in )
        update!( Ham, Rhoe_in )
        evals = diag_LOBPCG!( Ham, psi, Ham.precKin,
                              tol=ethr_last,
                              Nstates_conv=Ham.electrons.Nstates_occ )
        psi[:] = psi[:]/sqrt(dVol) # renormalize

        E_f, Ham.energies.mTS =
        update_Focc!( Ham.electrons.Focc, smear_fermi, smear_fermi_entropy,
                      evals, Float64(Ham.electrons.Nelectrons), kT )

        calc_rhoe!( Ham, psi, Rhoe_out )
        # Make sure no negative rhoe
        for ip in 1:length(Rhoe_out)
            if Rhoe_out[ip] < eps()
                println("Negative Rhoe is encountered")
                Rhoe_out[ip] = 0.0
            end
        end
        # Renormalize
        integ_rho = sum(Rhoe_out)*dVol
        #println("integ_rho (before renormalized) = ", integ_rho)
        for ip in 1:length(Rhoe_out)
            Rhoe_out[ip] = Ham.electrons.Nelectrons/integ_rho * Rhoe_out[ip]
        end
        #println("integ Rhoe before mix = ", sum(Rhoe_out)*dVol)
        
        calc_energies!( Ham, psi )
        Etot = sum(Ham.energies)
        dEtot = abs(Etot - Etot_old)
        Etot_old = Etot
        @printf("Etot: %18.10f %18.10e\n", Etot, dEtot)

        return betamix*Rhoe_out + (1 - betamix)*Rhoe_in
        #return Rhoe_out
    end

    NiterMax = 100
    println("\nNLsolve starting: \n")
    fpres = NLsolve.nlsolve(x -> density_map(x) - x, Rhoe,
            m=5, method=:anderson, xtol=1e-6, ftol=0.0, show_trace=true,
            iterations=NiterMax )
    println("\nNLsolve is converged in iter: ", fpres.iterations)

    println()
    calc_energies!( Ham, psi )
    println(Ham.energies)

    return
end