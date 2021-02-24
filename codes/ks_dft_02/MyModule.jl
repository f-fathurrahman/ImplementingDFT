module MyModule

using Printf
using LinearAlgebra
using SparseArrays
using FFTW

using AlgebraicMultigrid

include("../common/constants.jl")

export Ry2eV, ANG2BOHR

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

include("../common/GVectors.jl")
export GVectors

include("../common/ILU0Preconditioner.jl")

include("../common/diag_Emin_PCG.jl")
include("../common/diag_davidson.jl")
include("../common/diag_LOBPCG.jl")
export diag_LOBPCG!, diag_davidson!, diag_Emin_PCG!

include("../common/ortho_sqrt.jl")
export ortho_sqrt, ortho_sqrt!

include("../common/supporting_functions.jl")
export speye, meshgrid

include("../poisson_3d/Poisson_solve_PCG.jl")
export Poisson_solve_PCG

include("../poisson_3d/Poisson_solve_DAGE.jl")
export PoissonSolverDAGE, Poisson_solve_DAGE

include("../poisson_3d/Poisson_solve_fft.jl")
export PoissonSolverFFT, Poisson_solve_fft

export Poisson_solve

include("../common/Ylm_real.jl")

include("../common/Atoms.jl")
export Atoms

include("../common/calc_strfact.jl")
export calc_strfact, calc_strfact_shifted

include("../common/calc_dr_periodic.jl")
export calc_dr_periodic

include("Energies.jl")
export Energies

include("../common/PsPot_GTH.jl")
export PsPot_GTH, PsPot_GTH_octopus, eval_Vloc_R, eval_proj_R

#include("PsPotNL.jl")
include("PsPotNLSparse.jl")
export PsPotNL, check_betaNL_norm

include("Electrons.jl")
export Electrons

include("Hamiltonian.jl")
export Hamiltonian, update!
export op_H, op_V_Ps_nloc, op_K, op_V_Ps_loc, op_V_loc
export op_V_Ps_nloc!

include("calc_rhoe.jl")
export calc_rhoe, calc_rhoe!

include("../common/LDA_VWN.jl")
export excVWN, excpVWN

include("../common/calc_E_NN_ewald.jl")
include("calc_energies.jl")
export calc_E_NN, calc_energies!

end
