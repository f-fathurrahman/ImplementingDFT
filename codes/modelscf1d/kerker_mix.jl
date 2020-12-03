function  kerker_mix(res, mixOpts::KerkerMixOptions)
    #The Kerker mixing, following the parameter setup of elliptic
    #preconditioner.
    #
    #Lin Lin
    #Last revision: 5/6/2012
    return kerker_mix(res, mixOpts.KerkerB, mixOpts.kmult,
                      mixOpts.YukawaK, mixOpts.epsil0)

end


function  kerker_mix(res, KerkerB::Float64, kmult, YukawaK::Float64, epsil0::Float64)
    #The Kerker mixing, following the parameter setup of elliptic
    #preconditioner.
    #
    #Lin Lin
    #Last revision: 5/6/2012
    resfft = fft(res)
    resfft = ( epsil0 / (4*pi) * (kmult .+ YukawaK^2) )./
        ( KerkerB .+ epsil0 / (4*pi) * (kmult .+ YukawaK^2) ) .* resfft
    resprec = real(ifft(resfft))
    resprec = resprec .- mean(resprec)
    return resprec

end
