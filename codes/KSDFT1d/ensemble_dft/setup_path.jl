push!(LOAD_PATH, "../")

using Revise, Infiltrator
import Random
using Printf
using LinearAlgebra
using Serialization
using KSDFT1d

includet("../utilities.jl")
includet("../BroydenMixer.jl")
