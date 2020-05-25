include("../FD3d/FD3dGrid.jl")
include("../FD3d/build_nabla2_matrix.jl")

include("../LF3d/LF3dGrid.jl")
include("../LF3d/build_nabla2_matrix.jl")

include("../common/supporting_functions.jl")

include("Poisson_solve_CG.jl")
include("Poisson_solve_PCG.jl")