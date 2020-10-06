mutable struct LinearMixer
    betamix::Float64
end

function do_mix!(mixer::LinearMixer, x, x_new)
    mix_linear!(x, x_new, mixer.betamix)
    return
end

function mix_linear!(x, x_new, betamix::Float64)
    x = betamix*x_new + (1-betamix)*x
    return
end