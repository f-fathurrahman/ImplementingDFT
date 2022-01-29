module KSDFT1d

using Printf
using LinearAlgebra

include("INC_sch_1d.jl")

include("FD1dGrid.jl")
export FD1dGrid

include("Atoms1d.jl")
export Atoms1d

include("Electrons.jl")
export Electrons

include("smearing.jl")
export smear_fermi, smear_fermi_entropy

include("occupations.jl")
export update_Focc!

include("XCCalculator.jl")
include("Libxc_old.jl")
include("LDA_libxc_1d.jl")

include("Hamiltonian1d.jl")
export Hamiltonian1d

include("Poisson_solve_sum.jl")
export Poisson_solve_sum!

include("calc_rhoe.jl")
export calc_rhoe!

include("calc_energies.jl")
export calc_Ekin, calc_Enn

end # module