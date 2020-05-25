include("../FD1d/build_D2_matrix_3pt.jl")
include("../FD1d/build_D2_matrix_5pt.jl")
include("../FD1d/build_D2_matrix_7pt.jl")
include("../FD1d/build_D2_matrix_9pt.jl")
include("../FD1d/build_D2_matrix_11pt.jl")

include("../FD1d/build_D2_matrix_p_3pt.jl")
include("../FD1d/build_D2_matrix_p_5pt.jl")
include("../FD1d/build_D2_matrix_p_7pt.jl")
include("../FD1d/build_D2_matrix_p_9pt.jl")
include("../FD1d/build_D2_matrix_p_11pt.jl")

const ⊗ = kron

# XXX: Use metaprogramming?
function decide_D2_func1d( order, pbc::Bool )
    if order == 3
        if pbc
            return build_D2_matrix_p_3pt
        else
            return build_D2_matrix_3pt
        end
    
    elseif order == 5
        if pbc
            return build_D2_matrix_p_5pt
        else
            return build_D2_matrix_5pt
        end

    elseif order == 7
        if pbc
            return build_D2_matrix_p_7pt
        else
            return build_D2_matrix_7pt
        end

    elseif order == 11
        if pbc
            return build_D2_matrix_p_11pt
        else
            return build_D2_matrix_11pt
        end

    else # Default order=9
        if pbc
            return build_D2_matrix_p_9pt
        else
            return build_D2_matrix_9pt
        end

    end

end

function build_nabla2_matrix( grid::FD3dGrid; stencil_order=3 )

    func_1d = decide_D2_func1d( stencil_order, grid.pbc[1] )
    D2x = func_1d(grid.Nx, grid.hx)

    func_1d = decide_D2_func1d( stencil_order, grid.pbc[2] )    
    D2y = func_1d(grid.Ny, grid.hy)

    func_1d = decide_D2_func1d( stencil_order, grid.pbc[3] )
    D2z = func_1d(grid.Nz, grid.hz)

    IIx = speye(grid.Nx)
    IIy = speye(grid.Ny)
    IIz = speye(grid.Nz)

    ∇² = D2x⊗IIy⊗IIz + IIx⊗D2y⊗IIz + IIx⊗IIy⊗D2z 

    return ∇²

end