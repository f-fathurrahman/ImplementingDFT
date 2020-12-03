module ModelSCF1d

include("Atoms.jl")
export Atoms

include("scfOptions.jl")
export EigensolverOptions, SCFOptions

include("anderson_mix.jl")
export AndersonMixOptions

include("kerker_mix.jl")

include("pseudocharge.jl")
include("Hamiltonian.jl")
include("hartree_pot_bc.jl")

export Hamiltonian, init_pot!, scf_potmix!

include("getocc.jl")


#include("Integrators.jl")

end
