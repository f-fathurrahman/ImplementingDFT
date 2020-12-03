function lobpcg_sep(H, X0, prec, k::Int64;
                    tol::Float64=1e-8,
                    maxiter::Int64=50,
                    verbose::Bool=false)
# [X, lambda, iter] = lobpcg_sep(H, X0, , maxiter, X0, k, flag)
#
# Input:
#   - funA and prec are function handles for A and T, respectively.
#     They both need to accept multiple right hand sides.
#   - prec should take two additional arguments:
#     current Ritz values and the corresponding residual norms
#   - tol is the relative tolerance.
#   - X0 is the initial guess, which is assumed to be orthonormal
#     (this is an Wavefun object ).
#   - k is the number of desired eigenpairs.
#     k cannot exceed the number of columns of X.
#   - flag specifies whether debugging information is printed (1) or not (0)
#
# Output:
#   - X is a matrix that contains the normalized eigenvectors.
#   - lambda is a vector that contains the eigenvalues.
#   - iter is the number of iterations

X = X0;
n = size(X)[1];
tol_deflation = min(tol, 1024*eps(1.0));
AX = H*X;
Ap = X'*AX;

# symmetrizing Ap
Ap = (Ap' + Ap)/2; ##
(lambda, V) = eigen(Ap);
X = X*V;
AX = AX*V;
P = [];
AP = [];
D = diagm(0 => lambda);

# absolute residual
absr = zeros(k);
# relative residual
relr = zeros(k);

iteration = 1;
for iter = 1:maxiter
    # Check convergence
    R = AX - X*D;
    # number of converged eigenvectors
    nconv = 0;
    for j = 1:k
        absr[j] = norm(R[:, j], 2);
        relr[j] = absr[j]/(norm(AX[:, j], 2) + abs(lambda[j]));
        # Deflation
        if nconv == j-1 && relr[j] <= tol_deflation
            nconv = j;
        end
    end
    maxr = maximum(relr[1:k]);
    if verbose
        println("iter = ", iter);
        for j = 1:k
            println("  theta(",j,") = " ,lambda[j]," ,  res(",j,") = ",relr[j]);
        end
    end
    if maxr <= tol
        break;
    end
    # Apply the preconditioner
    W = zeros(size(R)[1], k-nconv);
    for ii = k:-1:nconv+1
        W[:,ii-nconv] = prec(R[:, ii]);
    end
    # Orthogonalization
    if ~isempty(P)
        # making sure that we are not concatenating different
        U = hcat(X, P);
    else
        U = X;
    end
    # Projecting the preconditioned residual
    W = W - U*(U'*W);
    W = qr(W).Q[:,1:size(W)[2]];
    AW = H*W;
    # Recalculate AX and AP in every 20 iterations
    if mod(iter, 20) == 19
        AX = H*X; AP = H*P;
    end
    # Projection
    U = hcat(U, W);
    if ~isempty(AP)
        # taking in acount if AP hasn't been initialized as a
        # wavefun yet
        AU = hcat(AX, AP, AW);
    else
        AU = hcat(AX, AW);
    end
    Ap = U'*AU; Ap = (Ap' + Ap)/2; ##
    (lambda, V) = eigen(Ap);
    # Update the basis
    Q = HL_orth(k, V);
    X  = U*Q[:, 1:k];
    P  = U*Q[:, k+1:2*k];
    AX = AU*Q[:, 1:k];
    AP = AU*Q[:, k+1:2*k];
    lambda = lambda[1:k];
    D = diagm(0 =>lambda)

    iteration = iter
end

X = X[:, 1:k];

return (lambda, X, iteration)

end

# qwe define a parametric function action on 2D arrays of any type
function HL_orth(n::Int64, Z::Array{T, 2}) where T <: Number
# Q = HL_orth(n, Z)
#
# Input:
#   - Z is an mxm unitary matrix.
#   - n is an integer satisfying 0 < 2n <= m.
#
# Output:
#   - Q is an mx2n matrix satisfying Q'*Q = I.
#
# Remark:
#   Z is partitioned as Z = [Z1; Z2], where Z1 is nxn.
#   Span(Q) should contain Span([Z1, 0; Z2, Z2]).
#

    m  = size(Z)[1];
    Q0 = qr(Z[1:n, n+1:m]').Q[:,1:m-n];
    Q  = hcat(Z[1:m, 1:n], Z[1:m, n+1:m]*Q0);

    return Q
end
