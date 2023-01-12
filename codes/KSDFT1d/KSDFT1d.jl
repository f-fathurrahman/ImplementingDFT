module KSDFT1d

using Printf
using LinearAlgebra

include("../common/constants.jl")

export Ry2eV, Ha2eV, eV2Ha, ANG2BOHR, BOHR2ANG

include("INC_sch_1d.jl")

include("FD1dGrid.jl")
export FD1dGrid

include("Atoms1d.jl")
export Atoms1d

include("Electrons.jl")
export Electrons

include("smearing.jl")
export smear_fermi, smear_fermi_entropy, smear_fermi_prime, grad_smear

include("occupations.jl")
export update_Focc!

include("XCCalculator.jl")
export XCCalculator

include("Energies.jl")
export Energies

include("Hamiltonian1d.jl")
export Hamiltonian1d
export op_H

include("Poisson_solve_sum.jl")
export Poisson_solve_sum!

include("calc_rhoe.jl")
export calc_rhoe!

include("calc_energies.jl")
export calc_E_kin, calc_E_NN

include("Libxc_old.jl")
include("LDA_libxc_1d.jl")
export calc_epsxc_1d, calc_Vxc_1d

end # module