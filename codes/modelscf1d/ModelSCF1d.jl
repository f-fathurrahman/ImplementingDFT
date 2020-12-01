module ModelSCF1d

#export Hamiltonian, Atoms, andersonMixOptions
#export eigOptions, scfOptions, init_pot!, scf!, get_force!
#export hartree_pot_bc, time_evolution, velocity_verlet, update_pot!

include("Atoms.jl")
export Atoms

#include("scfOptions.jl")
#include("anderson_mix.jl")
#include("kerker_mix.jl")
#include("Ham.jl")
#include("hartree_pot_bc.jl")
#include("pseudocharge.jl")
#include("getocc.jl")
#include("Integrators.jl")

end
