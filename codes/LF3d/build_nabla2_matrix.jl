include("../LF1d/build_D2_matrix_LF1d_c.jl")
include("../LF1d/build_D2_matrix_LF1d_sinc.jl")

const ⊗ = kron

function build_nabla2_matrix( grid::LF3dGrid )

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

    if grid.types[3] == :sinc
        D2z = build_D2_matrix_LF1d_sinc( grid.z, grid.hz, grid.Nz )
    else
        D2z = build_D2_matrix_LF1d_c( grid.Lz, grid.Nz )
    end

    IIx = speye(grid.Nx)
    IIy = speye(grid.Ny)
    IIz = speye(grid.Nz)

    ∇2 = D2x⊗IIy⊗IIz + IIx⊗D2y⊗IIz + IIx⊗IIy⊗D2z 

    return ∇2

end