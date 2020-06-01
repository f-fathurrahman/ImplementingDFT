function diag_LOBPCG!( Ham, X, prec;
                       tol=1e-5, NiterMax=100, verbose=false,
                       verbose_last=false, Nstates_conv=0 )

    Nbasis = size(X)[1]
    Nstates = size(X)[2]

    if Nstates_conv == 0
        Nstates_conv = Nstates
    end

    HX = Ham*X

    nconv = 0
    iter = 1
    resnrm = ones(Nstates)

    P = zeros(Float64,Nbasis,Nstates)
    HP = zeros(Float64,Nbasis,Nstates)

    evals = zeros(Float64,Nstates)
    devals = ones(Float64,Nstates)
    evals_old = copy(evals)

    SMALL = 10.0*eps()

    IS_CONVERGED = false

    for iter = 1:NiterMax

        # Rayleigh quotient
        XHX = Symmetric(X'*HX)
        evals = eigvals(XHX)

        devals = abs.( evals - evals_old )
        evals_old = copy(evals)
        
        # Residuals
        R = HX - X*XHX
        for j = 1:Nstates
            resnrm[j] = norm( @view R[:,j] )
        end

        nconv = length( findall( devals .< tol ) )
        
        ilock = findall( devals .<= SMALL )
        iact  = findall( devals .> SMALL )

        Nlock = length(ilock)
        Nact = length(iact)

        if verbose
            @printf("LOBPCG iter = %8d, nconv=%d\n", iter, nconv)
            for ist = 1:Nstates
                @printf("evals[%3d] = %18.10f, devals = %18.10e\n", ist, evals[ist], devals[ist] )
            end
            println("Nlock = ", Nlock, " Nact = ", Nact)
        end

        # Check for convergence
        if nconv >= Nstates_conv
            IS_CONVERGED = true
            if verbose || verbose_last
                @printf("LOBPCG convergence: nconv = %d in %d iterations\n", nconv, iter)
            end
            break
        end

        # apply preconditioner
        
        W = copy(R)
        for i in 1:size(W,2)
            @views ldiv!(prec, W[:,i])
        end

        if Nlock > 0
            Xlock = @view X[:,ilock]

            #Ylock = Kprec( Ham.ik, pw, Xlock )
            Ylock = copy(Xlock)
            for i in 1:size(Ylock,2)
                @views ldiv!(prec, Ylock[:,i])
            end

            Slock = Xlock'*Xlock
            for j = 1:Nact
                #W[:,iact[j]] = W[:,iact[j]] - Ylock*( Slock\( Xlock'*W[:,iact[j]] ) )
                ww = @view W[:,iact[j]]
                ww = ww - Ylock*( Slock\( Xlock'*ww ) )
            end
            #
            X[:,iact] = X[:,iact] - Xlock*(Xlock'*X[:,iact])
            W[:,iact] = W[:,iact] - Xlock*(Xlock'*W[:,iact])
            #
            HX[:,iact] = Ham*X[:,iact]
            HW = Ham*W[:,iact]
            #
            HX[:,iact] = HX[:,iact] - Xlock*( Xlock'*HX[:,iact] )
            HW = HW - Xlock*(Xlock'*HW)
            #
            C = W[:,iact]'*W[:,iact] #C = (C+C')/2;
            R = (cholesky(Symmetric(C))).U
            W[:,iact] = W[:,iact]/R
            HW = HW/R
            #
            Q  = [X[:,iact] W[:,iact]]
            HQ = [HX[:,iact] HW]

            T = Symmetric(Q'*(HQ))
            G = Symmetric(Q'*Q)
            sd, S = eigen( T, G ) # evals, evecs

            U = S[:,1:Nact]
            Xact = Q*U;
            X[:,:] = [Xlock Xact]
            T = Symmetric( X'* (Ham*X) ) #; T = (T+T')/2;
            evalT, evecsT = eigen(T);

            X[:,:] = X*evecsT
            HX = Ham*X

        else
            # nlock == 0
            HW = Ham*W
            #
            C  = W'*W
            #C = 0.5*( C + C' )
            #
            R  = (cholesky(Symmetric(C))).U
            W  = W/R
            HW = HW/R
            #
            Q  = [X W]
            HQ = [HX HW]
            if iter > 1
                Q  = [Q P]
                HQ = [HQ HP]
            end

            T = Symmetric(Q'*(HQ))
            G = Symmetric(Q'*Q)
            sd, S = eigen(T, G) # evals, evecs
            #
            U = S[:,1:Nstates]
            X[:,:] = Q*U
            HX = HQ*U
            if iter > 1
                set2 = Nstates+1:2*Nstates
                set3 = 2*Nstates+1:3*Nstates
                P  = W*U[set2,:] + P*U[set3,:]
                HP = HW*U[set2,:] + HP*U[set3,:]
                C = P'*P
                #C = 0.5*(C + C')
                R = (cholesky(Symmetric(C))).U
                P = P/R
                HP = HP/R
            else
                P  = copy(W)
                HP = copy(HW)
            end

        end  # if nlock 

        iter = iter + 1

        verbose && flush(stdout)

    end

    if !IS_CONVERGED && verbose
        @printf("\nWARNING: LOBPCG is not converged after %d iterations\n", NiterMax)
    end

    S = X'*HX
    S = Symmetric(0.5*(S+S'))
    evals, Q = eigen(S)
    X[:,:] = X[:,:]*Q
    if verbose_last || verbose
        @printf("\nEigenvalues from diag_LOBPCG:\n\n")
        for ist = 1:Nstates
            @printf("evals[%3d] = %18.10f, devals = %18.10e\n", ist, evals[ist], devals[ist] )
        end
    end
    return evals
end
