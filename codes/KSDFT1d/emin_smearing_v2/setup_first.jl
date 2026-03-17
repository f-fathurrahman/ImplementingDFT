push!(LOAD_PATH, "../")

using Revise
using Printf
using LinearAlgebra
using KSDFT1d
import Random
using Serialization

includet("../utilities.jl")
includet("../scf/solve_scf_01.jl")
