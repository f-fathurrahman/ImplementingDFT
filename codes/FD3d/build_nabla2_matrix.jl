include("../FD1d/build_D2_matrix_3pt.jl")
include("../FD1d/build_D2_matrix_5pt.jl")
include("../FD1d/build_D2_matrix_7pt.jl")
include("../FD1d/build_D2_matrix_9pt.jl")
include("../FD1d/build_D2_matrix_11pt.jl")

const ⊗ = kron

function build_nabla2_matrix( fdgrid::FD3dGrid; func_1d=build_D2_matrix_3pt )

    D2x = func_1d(fdgrid.Nx, fdgrid.hx)
    D2y = func_1d(fdgrid.Ny, fdgrid.hy)
    D2z = func_1d(fdgrid.Nz, fdgrid.hz)

    IIx = speye(fdgrid.Nx)
    IIy = speye(fdgrid.Ny)
    IIz = speye(fdgrid.Nz)

    ∇² = D2x⊗IIy⊗IIz + IIx⊗D2y⊗IIz + IIx⊗IIy⊗D2z 

    return ∇²

end