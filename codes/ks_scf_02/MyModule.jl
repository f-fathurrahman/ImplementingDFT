module MyModule

using Printf
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using IncompleteLU
using AlgebraicMultigrid

const Ry2eV = 13.605693009  # Ry to eV
const ANG2BOHR = 1.8897261254578281  # angstrom to bohr

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

include("../common/Ylm_real.jl")

include("../common/Atoms.jl")
include("../common/Atoms_io.jl")
include("../common/Atoms_utils.jl")
export Atoms

include("Energies.jl")
export Energies

include("PsPot_GTH.jl")
export PsPot_GTH, eval_Vloc_R, eval_proj_R

include("PsPotNL.jl")
export PsPotNL, calc_betaNL_psi

include("Electrons.jl")
export Electrons

include("Hamiltonian.jl")
export Hamiltonian, update!

include("calc_rhoe.jl")
export calc_rhoe, calc_rhoe!

include("../common/LDA_VWN.jl")
export excVWN, excpVWN

include("calc_energies.jl")
export calc_E_NN, calc_energies!

end
