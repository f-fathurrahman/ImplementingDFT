include("../LF1d/build_D2_matrix_LF1d_c.jl")
include("../LF1d/build_D2_matrix_LF1d_sinc.jl")

const ⊗ = kron

function build_nabla2_matrix( lfgrid::LF3dGrid )

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

    if lfgrid.type_z == :sinc
        D2z = build_D2_matrix_LF1d_sinc( lfgrid.z, lfgrid.hz, lfgrid.Nz )
    else
        D2z = build_D2_matrix_LF1d_c( lfgrid.Lz, lfgrid.Nz )
    end

    IIx = speye(lfgrid.Nx)
    IIy = speye(lfgrid.Ny)
    IIz = speye(lfgrid.Nz)

    ∇2 = D2x⊗IIy⊗IIz + IIx⊗D2y⊗IIz + IIx⊗IIy⊗D2z 

    return ∇2

end