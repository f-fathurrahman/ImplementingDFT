module MyModule

using Printf
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using IncompleteLU
using AlgebraicMultigrid

include("../FD2d/FD2dGrid.jl")
export FD2dGrid

include("../FD2d/build_nabla2_matrix.jl")
export build_nabla2_matrix
export build_D2_matrix_3pt,
       build_D2_matrix_5pt,
       build_D2_matrix_7pt,
       build_D2_matrix_9pt,
       build_D2_matrix_11pt

include("../LF2d/LF2dGrid.jl")
export LF2dGrid

include("../LF2d/build_nabla2_matrix.jl")
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

include("Poisson_solve_sum.jl")
export Poisson_solve_sum

include("Electrons.jl")
export Electrons

include("Energies.jl")
export Energies

include("Hamiltonian.jl")
export Hamiltonian, update!, op_H

include("calc_rhoe.jl")
export calc_rhoe, calc_rhoe!

include("LDA_2d.jl")
export calc_E_xc_2d, calc_V_xc_2d, calc_E_xc_V_xc_2d

include("calc_energies.jl")
export calc_energies!

end
