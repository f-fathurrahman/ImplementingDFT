struct EigensolverOptions
    # structure to encapsulate all the different options for the scf iterations
    # the main idea is to define a Hamiltonian, and then run scf on it; however,
    # the Hamiltoinian shouldn't have the scf options embeded
    eigstol::Float64
    eigsiter::Int64
    eigmethod::AbstractString
end

function EigensolverOptions()
    return EigensolverOptions(1e-8, 100, "eigs")
end


# defining abstract type for all the different mixing options
abstract type MixingOptions
end

mutable struct SimpleMixOptions <: MixingOptions
    betamix::Float64
end

# default constructor
function SimpleMixOptions()
    return SimpleMixOptions(0.5)
end

mutable struct AndersonMixOptions <: MixingOptions
    ymat::Array{Float64,2}
    smat::Array{Float64,2}
    betamix::Float64
    mixdim::Int64
    iter::Int64
end

function AndersonMixOptions(Ns, betamix, mixdim)
    ymat = zeros(Ns, mixdim)
    smat = zeros(Ns, mixdim)
    return AndersonMixOptions(ymat, smat, betamix[1], mixdim, 1)
end

mutable struct AndersonPrecMixOptions <: MixingOptions
    ymat::Array{Float64,2}
    smat::Array{Float64,2}
    betamix::Float64
    mixdim::Int64
    iter::Int64
    prec::Function
    precArgs
end


function AndersonPrecMixOptions(Ns, betamix, mixdim, prec, precArgs)
    ymat = zeros(Ns, scfOpts.mixdim);
    smat = zeros(Ns, scfOpts.mixdim);
    new(ymat,smat,betamix[1],mixdim,1, prec, precArgs)
end


mutable struct KerkerMixOptions <: MixingOptions
    betamix::Float64
    KerkerB::Float64
    kx::Array{Float64,2}
    YukawaK::Float64
    epsil0::Float64
end

function updateMix!( mixOpts::AndersonMixOptions, ii )
    mixOpts.iter = ii
    return
end

struct SCFOptions
    SCFtol::Float64
    scfiter::Int64
    eigOpts::EigensolverOptions
    mixOpts::MixingOptions
end

# TODO: this Initializing is ugly, we need to use kwargs to properly initialize
function SCFOptions()
    eigOpts = EigensolverOptions()
    mixOpts = SimpleMixOptions()
    return SCFOptions(1e-7,100, eigOpts, mixOpts);
end
    
function SCFOptions(eigOpts::EigensolverOptions, mixOpts::MixingOptions)
    return SCFOptions(1e-7,100,  eigOpts, mixOpts)
end
