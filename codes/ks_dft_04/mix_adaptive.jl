mutable struct AdaptiveLinearMixer
    betamix::Float64
    betamax::Float64
    v::Array{Float64,1}
    df::Array{Float64,1}
end

# Constructor
function AdaptiveLinearMixer(
    Npoints::Int64,
    betamix::Float64, betamax::Float64;
    Nspin=1
)
    v = betamix*ones(Float64, Npoints*Nspin)
    df = zeros(Float64, Npoints*Nspin)
    return AdaptiveLinearMixer(betamix, betamax, v, df)
end

# Result is put in v
function do_mix!(mixer::AdaptiveLinearMixer, v, v_new)
    mix_adaptive!( v, v_new,
        mixer.betamix, mixer.v, mixer.df, betamax=mixer.betamax
    )
    return
end

function mix_adaptive!( mu, nu, beta0::Float64, beta, f; betamax=0.8 )

    Npts = length(mu)
    
    for i = 1:Npts
        
        t = nu[i] - mu[i]

        if t*f[i] >= 0.0
            beta[i] = beta[i] + beta0
            if beta[i] > betamax
                beta[i] = betamax
            end
        else
            beta[i] = 0.5*( beta[i] + beta0 )
        end
        
        f[i] = t
        
        mu[i] = beta[i]*nu[i] + ( 1.0 - beta[i] )*mu[i]

    end

    return

end