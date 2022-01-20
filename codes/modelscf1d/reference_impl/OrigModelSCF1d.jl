module OrigModelSCF1d

include("Atoms.jl")
export Atoms

include("scfOptions.jl")
export EigensolverOptions, SCFOptions
export AndersonMixOptions, SimpleMixOptions

include("anderson_mix.jl")
include("kerker_mix.jl")

include("pseudocharge.jl")
export pseudocharge

include("Hamiltonian.jl")
include("hartree_pot_bc.jl")
export hartree_pot_bc
export Hamiltonian, init_pot!, scf_potmix!, get_force!
export create_Hamiltonian

include("getocc.jl")


#include("Integrators.jl")

end
