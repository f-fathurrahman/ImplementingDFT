module MyModule

using Printf
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using IncompleteLU
using AlgebraicMultigrid

include("../../3d/FD3dGrid.jl")
export FD3dGrid

include("../../3d/build_nabla2_matrix.jl")
export build_nabla2_matrix
export build_D2_matrix_3pt,
       build_D2_matrix_5pt,
       build_D2_matrix_7pt,
       build_D2_matrix_9pt

include("../../diag_Emin_PCG.jl")
include("../../diag_davidson.jl")
include("../../diag_LOBPCG.jl")
export diag_LOBPCG!, diag_davidson!, diag_Emin_PCG!

include("../../ortho_sqrt.jl")
export ortho_sqrt, ortho_sqrt!

include("../../supporting_functions.jl")
export speye, meshgrid

include("../../3d_poisson/Poisson_solve_PCG.jl")
export Poisson_solve_PCG

include("Electrons.jl")
export Electrons

include("Energies.jl")
export Energies

include("Hamiltonian.jl")
export Hamiltonian, update!

include("calc_rhoe.jl")
export calc_rhoe, calc_rhoe!

include("../LDA_VWN.jl")
export excVWN, excpVWN

include("calc_energies.jl")
export calc_energies!

end
