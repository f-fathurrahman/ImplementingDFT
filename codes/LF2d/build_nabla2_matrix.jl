include("../LF1d/build_D2_matrix_LF1d_c.jl")
include("../LF1d/build_D2_matrix_LF1d_sinc.jl")

function build_nabla2_matrix( lfgrid::LF2dGrid )

    if lfgrid.type_x == :sinc
        D2x = build_D2_matrix_LF1d_sinc( lfgrid.x, lfgrid.hx, lfgrid.Nx )
    else
        D2x = build_D2_matrix_LF1d_c( lfgrid.Lx, lfgrid.Nx )
    end

    if lfgrid.type_y == :sinc
        D2y = build_D2_matrix_LF1d_sinc( lfgrid.y, lfgrid.hy, lfgrid.Ny )
    else
        D2y = build_D2_matrix_LF1d_c( lfgrid.Ly, lfgrid.Ny )
    end

    ∇2 = kron(D2x, speye(lfgrid.Ny)) + kron(speye(lfgrid.Nx), D2y)
    return ∇2

end