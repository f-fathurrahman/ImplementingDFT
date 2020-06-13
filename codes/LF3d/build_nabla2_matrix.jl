include("../LF1d/build_D2_matrix_LF1d_c.jl")
include("../LF1d/build_D2_matrix_LF1d_sinc.jl")
include("../LF1d/build_D2_matrix_LF1d_p.jl")

const ⊗ = kron

function build_nabla2_matrix( grid::LF3dGrid )

    if grid.types[1] == :sinc
        D2x = build_D2_matrix_LF1d_sinc( grid.x, grid.hx, grid.Nx )
    elseif grid.types[1] == :C
        D2x = build_D2_matrix_LF1d_c( grid.Lx, grid.Nx )
    elseif grid.types[1] == :P
        D2x = build_D2_matrix_LF1d_p( grid.Lx, grid.Nx )
    else
        error("Unsupported types[1] = ", types[1])
    end


    if grid.types[2] == :sinc
        D2y = build_D2_matrix_LF1d_sinc( grid.y, grid.hy, grid.Ny )
    elseif grid.types[2] == :C
        D2y = build_D2_matrix_LF1d_c( grid.Ly, grid.Ny )
    elseif grid.types[2] == :P
        D2y = build_D2_matrix_LF1d_p( grid.Ly, grid.Ny )
    else
        error("Unsupported types[2] = ", types[2])
    end

    if grid.types[3] == :sinc
        D2z = build_D2_matrix_LF1d_sinc( grid.z, grid.hz, grid.Nz )
    elseif grid.types[3] == :C
        D2z = build_D2_matrix_LF1d_c( grid.Lz, grid.Nz )
    elseif grid.types[3] == :P
        D2z = build_D2_matrix_LF1d_p( grid.Lz, grid.Nz )
    else
        error("Unsupported types[3] = ", types[3])
    end

    IIx = speye(grid.Nx)
    IIy = speye(grid.Ny)
    IIz = speye(grid.Nz)

    ∇2 = D2x⊗IIy⊗IIz + IIx⊗D2y⊗IIz + IIx⊗IIy⊗D2z 

    return ∇2

end