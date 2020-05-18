include("../FD1d/build_D2_matrix_3pt.jl")
include("../FD1d/build_D2_matrix_5pt.jl")
include("../FD1d/build_D2_matrix_7pt.jl")
include("../FD1d/build_D2_matrix_9pt.jl")
include("../FD1d/build_D2_matrix_11pt.jl")

function build_nabla2_matrix( fdgrid::FD2dGrid; func_1d=build_D2_matrix_3pt )
    Nx = fdgrid.Nx
    hx = fdgrid.hx
    Ny = fdgrid.Ny
    hy = fdgrid.hy
    
    D2x = func_1d(Nx, hx)
    D2y = func_1d(Ny, hy)

    ∇² = kron(D2x, speye(Ny)) + kron(speye(Nx), D2y)
    return ∇²
end