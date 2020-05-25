include("../LF1d/build_D2_matrix_LF1d_c.jl")
include("../LF1d/build_D2_matrix_LF1d_sinc.jl")

function build_nabla2_matrix( grid::LF2dGrid )

    if grid.types[1] == :sinc
        D2x = build_D2_matrix_LF1d_sinc( grid.x, grid.hx, grid.Nx )
    else
        D2x = build_D2_matrix_LF1d_c( grid.Lx, grid.Nx )
    end

    if grid.types[2] == :sinc
        D2y = build_D2_matrix_LF1d_sinc( grid.y, grid.hy, grid.Ny )
    else
        D2y = build_D2_matrix_LF1d_c( grid.Ly, grid.Ny )
    end

    ∇2 = kron(D2x, speye(grid.Ny)) + kron(speye(grid.Nx), D2y)
    return ∇2

end