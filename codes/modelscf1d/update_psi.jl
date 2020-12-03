function update_psi!(H::Hamiltonian, eigOpts::EigensolverOptions)
    # we need to add some options to the update
    # functio to solve the eigenvalue problem for a given rho and Vtot

    if eigOpts.eigmethod == "eigs"
        # TODO: make sure that eigs works with overloaded operators
        # TODO: take a look a the interface of eigs in Julia
        (ev,psi,nconv,niter,nmult,resid) = eigs(H,   # Hamiltonian
                                                nev=H.Neigs+4, # number of eigs
                                                which=:SR, # small real part
                                                ritzvec=true, # provide Ritz v
                                                tol=eigOpts.eigstol, # tolerance
                                                maxiter=eigOpts.eigsiter) #maxiter
        #assert(flag == 0)
        # TODO: this has a bug somewhere !!!! fix me!!!!
    elseif  eigOpts.eigmethod == "lobpcg_sep"
        # not working... to be fixed
        X0 = qr(rand(H.Ns, H.Neigs)).Q[:,1:H.Neigs]
        precond(x) = inv_lap(H,x)

        (ev,psi, iter) = lobpcg_sep(H, X0, precond, H.Neigs,
                                    tol=eigOpts.eigstol,
                                    maxiter=eigOpts.eigsiter)

    elseif  eigOpts.eigmethod == "eig"
        # we use a dense diagonalization
        A = create_Hamiltonian(H)
        # checkign that A is symetric
        @assert issymmetric(A)
        (ev, psi) = (eigen(A)...,)

    end

    # sorting the eigenvalues, eigs already providesd them within a vector
    ind = sortperm(ev)[1:H.Neigs]
    # updating the eigenvalues
    H.ev = ev[ind]
    # updating the eigenvectors
    H.psi = psi[:, ind]
end