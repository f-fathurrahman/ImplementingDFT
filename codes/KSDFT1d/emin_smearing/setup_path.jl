push!(LOAD_PATH, "../")

using Revise, Infiltrator
import Random
using Printf
using LinearAlgebra
using Serialization
using KSDFT1d

includet("system_defs_01.jl")
includet("system_defs_02.jl")
includet("Lfunc.jl")
includet("../utilities.jl")
include("gradients_psi_Haux.jl")