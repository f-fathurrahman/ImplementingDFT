using Arpack

# TODO: add a default setting for the scfOpts
function scf_potmix!(H::Hamiltonian, scfOpts::SCFOptions)

    VtoterrHist = zeros(scfOpts.scfiter)
    eigOpts = scfOpts.eigOpts
    mixOpts = scfOpts.mixOpts
    Nocc = round(Integer, sum(H.atoms.nocc) / H.nspin)

    for ii = 1:scfOpts.scfiter
        # solving the linear eigenvalues problem
        update_psi!(H, eigOpts)

        # update the electron density
        update_rho!(H, Nocc)

        # update the total potential, and compute the
        # differnce between the potentials
        Verr = update_vtot!(H, mixOpts)

        # save the error
        VtoterrHist[ii] = Verr

        @printf("iter = %3d Verr = %18.10f\n", ii, Verr)
        # test if the problem had already satiesfied the tolerance
        if scfOpts.SCFtol > Verr
            break
        end
    end

    return VtoterrHist[VtoterrHist.>0]
end