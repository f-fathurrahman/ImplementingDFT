module ModelSCF1d

using Printf
using LinearAlgebra
using FFTW

include("Atoms1d.jl")
export Atoms1d

include("GVectors1d.jl")
export GVectors1d

include("Electrons.jl")
export Electrons

include("smearing.jl")
export smear_fermi, smear_fermi_entropy

include("occupations.jl")
export update_Focc!

include("Hamiltonian1d.jl")
export Hamiltonian1d, get_matrix

include("PoissonYukawa1d_solve.jl")
export PoissonYukawa1d_solve!

end
