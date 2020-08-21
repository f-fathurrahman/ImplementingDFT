# used by diag_Emin_PCG
function calc_grad_evals!( Ham, psi, grad )
    Nbasis  = size(psi,1)
    Nstates = size(psi,2)
    #
    H_psi = Ham*psi
    for ist in 1:Nstates
        @views grad[:,ist] = H_psi[:,ist]
        for jst in 1:Nstates
            @views grad[:,ist] = grad[:,ist] - dot( psi[:,jst], H_psi[:,ist] ) * psi[:,jst]
        end
    end
    return grad
end


function diag_Emin_PCG!(
    Ham, X::Array{Float64,2}, prec;
    tol=1e-5,
    NiterMax=100,
    verbose=false,
    verbose_last=false,
    Nstates_conv=0,
    tol_ebands=1e-4,
    α_t=3e-5,
    conv_info=[0,0]
)

    Nbasis  = size(X,1)
    Nstates = size(X,2)

    if Nstates_conv == 0
        Nstates_conv = Nstates
    end

    #
    # Variabls for PCG
    #
    g = zeros(Float64,Nbasis,Nstates)
    d = copy(g)
    g_old = copy(g)
    d_old = copy(g)
    Kg = copy(g)
    Xc = copy(g)
    gt = copy(g)
    
    β = 0.0
    
    Ebands_old = 0.0

    Hr = Hermitian( X' * (Ham*X) )

    evals = eigvals(Hr)
    Ebands = sum(evals)
    
    evals_old = copy(evals)
    devals = ones(Nstates)

    IS_CONVERGED = false
    gKg_old = 1.0

    for iter = 1:NiterMax

        calc_grad_evals!( Ham, X, g )

        @views Kg[:,:] .= g[:,:]
        for i in 1:Nstates
            @views ldiv!(prec, Kg[:,i])
        end

        gKg = real(dot(g,Kg))
        if iter != 1
            gOldKg = real(dot(g_old,Kg))
            β = (gKg - gOldKg) / gKg_old
        end
        if β < 0.0
            β = 0.0
        end

        d = -Kg + β * d_old

        Xc = ortho_sqrt( X + α_t*d )
        calc_grad_evals!( Ham, Xc, gt )

        denum = real(dot(g-gt,d))
        if denum != 0.0
            α = abs( α_t*real(dot(g,d))/denum )
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
            conv_info[1] = nconv # no of converged eigevalues
            conv_info[2] = iter  # no of iterations to converge
            if verbose || verbose_last
                @printf("Convergence is achieved based on nconv\n")
            end
            break
        end

        d_old = copy(d)
        g_old = copy(g)
        gKg_old = gKg
        Ebands_old = Ebands

        verbose && flush(stdout)

    end

    if !IS_CONVERGED
        @printf("\nWARNING: diag_Emin_PCG is not converged after %d iterations\n", NiterMax)
    end

    ortho_sqrt!(X)
    Hr = Hermitian( X' * (Ham*X) )
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
