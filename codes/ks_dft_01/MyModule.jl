module MyModule

using Printf
using LinearAlgebra
using SparseArrays
using AlgebraicMultigrid

include("../FD3d/FD3dGrid.jl")
export FD3dGrid

include("../FD3d/build_nabla2_matrix.jl")
export build_nabla2_matrix
export build_D2_matrix_3pt,
       build_D2_matrix_5pt,
       build_D2_matrix_7pt,
       build_D2_matrix_9pt,
       build_D2_matrix_11pt

include("../LF3d/LF3dGrid.jl")
export LF3dGrid

include("../LF3d/build_nabla2_matrix.jl")
export build_D2_matrix_LF1d_c,
       build_D2_matrix_LF1d_sinc

include("../common/diag_Emin_PCG.jl")
include("../common/diag_davidson.jl")
include("../common/diag_LOBPCG.jl")
export diag_LOBPCG!, diag_davidson!, diag_Emin_PCG!

include("../common/ortho_sqrt.jl")
export ortho_sqrt, ortho_sqrt!

include("../common/supporting_functions.jl")
export speye, meshgrid

include("../common/constants.jl")
export Ry2eV, ANG2BOHR

include("../common/Atoms.jl")
export Atoms

include("../common/GVectors.jl")
export GVectors

include("../common/ILU0Preconditioner.jl")

include("../poisson_3d/Poisson_solve_PCG.jl")
include("../poisson_3d/Poisson_solve_fft.jl")
include("../poisson_3d/Poisson_solve_DAGE.jl")
export Poisson_solve_PCG, Poisson_solve_fft, Poisson_solve_DAGE
export Poisson_solve

include("Electrons.jl")
export Electrons

include("Energies.jl")
export Energies

include("Hamiltonian.jl")
export Hamiltonian, update!

include("calc_rhoe.jl")
export calc_rhoe, calc_rhoe!

include("../common/LDA_VWN.jl")
export excVWN, excpVWN

include("../common/XCCalculator.jl")
export XCCalculator

include("../common/XC_functionals_internal.jl")
include("../common/LDA_VWN_internal.jl")
export calc_epsxc_VWN, calc_Vxc_VWN

include("calc_energies.jl")
export calc_energies!

end
