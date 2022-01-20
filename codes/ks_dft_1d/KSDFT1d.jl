module KSDFT1d

using Printf
using LinearAlgebra

include("../FD1d/init_FD1d_grid.jl")
include("../FD1d/build_D2_matrix_3pt.jl")
include("../FD1d/build_D2_matrix_5pt.jl")
include("../FD1d/build_D2_matrix_7pt.jl")
include("../FD1d/build_D2_matrix_9pt.jl")
include("../FD1d/build_D2_matrix_11pt.jl")

include("../FD1d/init_FD1d_p_grid.jl")
#
include("../FD1d/build_D1_matrix_p_3pt.jl")
include("../FD1d/build_D1_matrix_p_5pt.jl")
include("../FD1d/build_D1_matrix_p_7pt.jl")
include("../FD1d/build_D1_matrix_p_9pt.jl")
include("../FD1d/build_D1_matrix_p_11pt.jl")
#
include("../FD1d/build_D2_matrix_p_3pt.jl")
include("../FD1d/build_D2_matrix_p_5pt.jl")
include("../FD1d/build_D2_matrix_p_7pt.jl")
include("../FD1d/build_D2_matrix_p_9pt.jl")
include("../FD1d/build_D2_matrix_p_11pt.jl")

include("FD1dGrid.jl")
export FD1dGrid

end # module