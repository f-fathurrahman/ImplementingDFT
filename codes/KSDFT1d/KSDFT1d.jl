module KSDFT1d

using Printf
using LinearAlgebra

# constants
#
# CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?rydhcev
const Ry2eV = 13.605693009  # Ry to eV
const Ha2eV = 2*Ry2eV
const eV2Ha = 1/Ha2eV

#
# CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
# 1/bohr
const ANG2BOHR = 1.8897261254578281  # angstrom to bohr
const BOHR2ANG = 1/ANG2BOHR

export Ry2eV, Ha2eV, eV2Ha, ANG2BOHR, BOHR2ANG

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
export XCCalculator

include("Hamiltonian1d.jl")
export Hamiltonian1d

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