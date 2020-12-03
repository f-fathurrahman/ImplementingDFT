function update_vtot!(H::Hamiltonian, mixOpts)
    # TODO here we need to implement the andersson mix
    # I added the signature

    (Vtotnew,Verr) = update_pot!(H)

    # TODO: add a swtich to use different kinds of mixing here
    betamix = mixOpts.betamix
    mixdim = mixOpts.mixdim
    ymat = mixOpts.ymat
    smat = mixOpts.smat
    iter = mixOpts.iter

    (Vtotmix,ymat,smat) = anderson_mix(H.Vtot, Vtotnew,
        betamix, ymat, smat, iter, mixdim)

    # they are already modified inside the function
    # mixOpts.ymat = ymat
    # mixOpts.smat = smat
    mixOpts.iter += 1

    # updating total potential
    H.Vtot = Vtotmix
    return Verr
end

function update_vtot!(H::Hamiltonian, mixOpts::KerkerMixOptions)
    # TODO here we need to implement the andersson mix
    # I added the signature

    (Vtotnew,Verr) = update_pot!(H)

    # println("using kerker mixing")
    # computign the residual
    res     = H.Vtot - Vtotnew

    # println("appliying the preconditioner")
    resprec = kerker_mix(res, mixOpts.KerkerB, mixOpts.kx,
                              mixOpts.YukawaK, mixOpts.epsil0)

    # println("performing the linear mixing ")
    Vtotmix =  H.Vtot - mixOpts.betamix * resprec

    # updating total potential
    # println("saving the total potentail ")
    H.Vtot = Vtotmix

    # println("returning the error")
    return Verr
end

function update_vtot!(H::Hamiltonian, mixOpts::AndersonPrecMixOptions)
    # TODO here we need to implement the andersson mix
    # I added the signature

    (Vtotnew,Verr) = update_pot!(H)

    betamix = mixOpts.betamix
    mixdim = mixOpts.mixdim
    ymat = mixOpts.ymat
    smat = mixOpts.smat
    iter = mixOpts.iter

    # TODO change this so this is properly modified inside the function
    (Vtotmix,ymat,smat) = prec_anderson_mix(H.Vtot,Vtotnew,
        betamix, ymat, smat, iter, mixdim, mixOpts.prec, mixOpts.precargs)

    # they are already modified inside the function
    # just to be sure we can leave them here
    mixOpts.ymat = ymat
    mixOpts.smat = smat
    mixOpts.iter += 1

    # updating total potential
    H.Vtot = Vtotmix
    return Verr
end