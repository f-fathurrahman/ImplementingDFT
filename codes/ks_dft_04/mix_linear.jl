mutable struct LinearMixer
    betamix::Float64
    iter::Int64 # needed for compatibility with other mixer
end

function LinearMixer(betamix)
    return LinearMixer(betamix, 0)
end

function do_mix!(mixer::LinearMixer, x, x_new)
    mix_linear!(x, x_new, mixer.betamix)
    return
end

function mix_linear!(x, x_new, betamix::Float64)
    Npts = length(x)
    for i in 1:Npts
        x[i] = betamix*x_new[i] + (1-betamix)*x[i]
    end
    return
end