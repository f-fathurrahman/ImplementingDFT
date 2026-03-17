push!(LOAD_PATH, "../")

using Revise
using Printf
using LinearAlgebra
using KSDFT1d
import Random

includet("../utilities.jl")

includet("system_defs_01.jl")
includet("common_emin_no_smearing.jl")
includet("solve_emin_CG_v01.jl")