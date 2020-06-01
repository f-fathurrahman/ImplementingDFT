function calc_grad_evals!( Ham, ψ, g, Hsub )
    Nstates = size(ψ,2)
    Hψ = Ham*ψ
    Hsub[:] = ψ' * Hψ
    g[:,:] = Hψ - ψ*Hsub
    return
end


function diag_Emin_PCG!( Ham, X, prec;
                         tol=1e-5,
                         NiterMax=100,
                         verbose=false,
                         verbose_last=false,
                         Nstates_conv=0,
                         tol_ebands=1e-4,
                         α_t=3e-5 )

    Nbasis  = size(X,1)
    Nstates = size(X,2)

    if Nstates_conv == 0
        Nstates_conv = Nstates
    end

    #
    # Variabls for PCG
    #
    Hsub = zeros(Float64,Nstates,Nstates)
    g = zeros(Float64,Nbasis,Nstates)
    d = zeros(Float64,Nbasis,Nstates)
    g_old = copy(d)
    d_old = copy(d)
    Kg = copy(d)
    Kg_old = copy(d)
    Xc = copy(d)
    gt = copy(d)
    
    β = 0.0
    
    Ebands_old = 0.0

    Hr = Symmetric( X' * (Ham*X) )

    evals = eigvals(Hr)
    Ebands = sum(evals)
    
    evals_old = copy(evals)
    devals = ones(Nstates)

    IS_CONVERGED = false

    for iter = 1:NiterMax

        calc_grad_evals!( Ham, X, g, Hsub )
        
        Kg[:] = g[:] # copy
        for i in 1:Nstates
            @views ldiv!(prec, Kg[:,i])
        end

        if iter != 1
            β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g_old).*Kg_old))
        end

        if β < 0.0
            β = 0.0
        end

        d = -Kg + β * d_old

        Xc = ortho_sqrt( X + α_t*d )
        calc_grad_evals!( Ham, Xc, gt, Hsub )

        denum = real(sum(conj(g-gt).*d))
        if denum != 0.0
            α = abs( α_t*real(sum(conj(g).*d))/denum )
        else
            α = 0.0
        end

        # Update wavefunction
        X[:,:] = X + α*d
        ortho_sqrt!(X)

        Hr = Symmetric( X' * ( Ham*X ) )
        evals = eigvals(Hr)
        Ebands = sum(evals)

        devals = abs.( evals - evals_old )
        evals_old = copy(evals)
        
        nconv = length( findall( devals .< tol ) )

        diffE = abs(Ebands-Ebands_old)

        if verbose
            @printf("CG step %8d = %18.10f %10.7e\n", iter, Ebands, diffE)
            for ist = 1:Nstates
                @printf("evals[%3d] = %18.10f, devals = %18.10e\n", ist, evals[ist], devals[ist] )
            end
            @printf("iter %d nconv = %d\n", iter, nconv)
        end
        if nconv >= Nstates_conv
            IS_CONVERGED = true
            if verbose || verbose_last
                @printf("Emin_PCG convergence: nconv = %d in %d iterations\n", nconv, iter)
            end
            break
        end

        g_old = copy(g)
        d_old = copy(d)
        Kg_old = copy(Kg)
        Ebands_old = Ebands

        if verbose flush(stdout) end

    end

    if !IS_CONVERGED
        @printf("\nWARNING: diag_Emin_PCG is not converged after %d iterations\n", NiterMax)
    end

    ortho_sqrt!(X)
    Hr = Symmetric( X' * (Ham*X) )
    evals, evecs = eigen(Hr)
    X[:,:] = X*evecs

    if verbose_last || verbose
        @printf("\nEigenvalues from diag_Emin_PCG:\n\n")
        for ist = 1:Nstates
            @printf("evals[%3d] = %18.10f devals = %18.10e\n", ist, evals[ist], devals[ist] )
        end
    end

    return evals
end
