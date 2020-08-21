"""
Assumption:
dot(X,X) = 1
"""
function diag_davidson!(
    Ham, X::Array{Float64,2}, prec;
    tol=1e-5, NiterMax=100, verbose=false,
    verbose_last=false, Nstates_conv=0,
    conv_info=[0,0]
)


    # get size info
    Nstates = size(X)[2]
    Nbasis  = size(X)[1]

    @assert(Nstates >= 1)

    if Nstates_conv == 0
        Nstates_conv = Nstates
    end

    evals    = zeros(Float64, Nstates)
    R        = zeros(Float64, Nbasis, Nstates)
    Hred     = zeros(Float64, 2*Nstates, 2*Nstates)
    Sred     = zeros(Float64, 2*Nstates, 2*Nstates)
    λ_red    = zeros(Float64, 2*Nstates)
    X_red    = zeros(Float64, 2*Nstates, 2*Nstates)
    res      = zeros(Float64, Nstates)
    res_norm = zeros(Float64, Nstates)

    devals    = zeros(Float64, Nstates)
    evals_old = zeros(Float64, Nstates)

    HX = Ham*X

    # Initial eigenvalues
    for ist = 1:Nstates
        evals[ist] = real( dot( X[:,ist], HX[:,ist] ) )
    end

    # Calculate residuals
    for ist = 1:Nstates
        for ig = 1:Nbasis
            R[ig,ist] = evals[ist]*X[ig,ist] - HX[ig,ist]
        end
        res[ist] = sqrt( dot( R[:,ist], R[:,ist] ) )
    end

    EPS = eps()

    set1 = 1:Nstates
    set2 = Nstates+1:2*Nstates

    evals_old = copy(evals)

    IS_CONVERGED = false

    for iter = 1:NiterMax

        res_norm[:] .= 1.0

        for ist = 1:Nstates
            if EPS < res[ist]
                res_norm[ist] = 1.0/res[ist]
            else
                println("res[ic] is too small")
            end
            for ig = 1:Nbasis
                R[ig,ist] = res_norm[ist] * R[ig,ist]
            end
        end

        #R = Kprec(Ham.ik, Ham.pw, R)
        for i in 1:Nstates
            @views ldiv!(prec, R[:,i])
        end

        HR = Ham*R

        # FIXME: Pull this outside the loop?
        if iter == 1
            Hred[set1,set1] = X' * HX
        else
            Hred[1:Nstates,1:Nstates] .= 0.0
            for ic = 1:Nstates
                Hred[ic,ic] = evals[ic] # use diagm ?
            end
        end

        Hred[set1,set2] = X' * HR
        Hred[set2,set2] = R' * HR
        Hred[set2,set1] = Hred[set1,set2]'

        Sred[set1,set1] = Matrix(Diagonal(ones(Float64,Nstates)))
        Sred[set1,set2] = X' * R
        Sred[set2,set2] = R' * R
        Sred[set2,set1] = Sred[set1,set2]'

        Hred = 0.5*(Hred + Hred')
        Sred = 0.5*(Sred + Sred')

        λ_red, X_red = eigen( Symmetric(Hred), Symmetric(Sred) )

        evals = λ_red[1:Nstates]

        devals = abs.( evals - evals_old )
        nconv = length( findall( devals .< tol ) )

        X[:,:]  = X  * X_red[set1,set1] + R  * X_red[set2,set1]
        HX = HX * X_red[set1,set1] + HR * X_red[set2,set1]

        # Calculate residuals
        for ist = 1:Nstates
            for ig = 1:Nbasis
                R[ig,ist] = evals[ist]*X[ig,ist] - HX[ig,ist]
            end
            res[ist] = sqrt( dot( R[:,ist], R[:,ist] ) )
        end

        if verbose
            @printf("\n")
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
                @printf("Davidson convergence: Nstates_conv in iter = %d\n", iter)
            end
            break
        end

        evals_old = copy(evals)
    end

    if !IS_CONVERGED
        @printf("\nWARNING: diag_davidson is not converged after %d iterations\n", NiterMax)
    end

    if verbose_last || verbose
        @printf("\nEigenvalues from diag_davidson:\n\n")
        for ist = 1:Nstates
            @printf("evals[%3d] = %18.10f, devals = %18.10e\n", ist, evals[ist], devals[ist] )
        end
    end

    return evals

end
