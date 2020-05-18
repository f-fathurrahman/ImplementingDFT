include("../LF1d/build_D2_matrix_LF1d_c.jl")

function build_nabla2_matrix( lfgrid::LF2dGrid )
    
    Lx = lfgrid.Lx
    Nx = lfgrid.Nx
    
    Ly = lfgrid.Ly
    Ny = lfgrid.Ny
    
    D2x = build_D2_matrix_LF1d_c(Lx, Nx)
    D2y = build_D2_matrix_LF1d_c(Lx, Ny)

    ∇2 = kron(D2x, speye(Ny)) + kron(speye(Nx), D2y)
    return ∇2
end