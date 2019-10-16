include("../1d/build_D2_matrix_3pt.jl")
include("../1d/build_D2_matrix_5pt.jl")
include("../1d/build_D2_matrix_7pt.jl")
include("../1d/build_D2_matrix_9pt.jl")

function build_nabla2_matrix( fdgrid::FD2dGrid; func_1d=build_D2_matrix_3pt )
    Nx = fdgrid.Nx
    hx = fdgrid.hx
    Ny = fdgrid.Ny
    hy = fdgrid.hy
    
    D2x = func_1d(Nx, hx)
    D2y = func_1d(Ny, hy)

    ∇2 = kron(D2x, speye(Ny)) + kron(speye(Nx), D2y)
    return ∇2
end