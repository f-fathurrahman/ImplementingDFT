speye(N::Int64) = sparse( Matrix(1.0I, N, N) )

"""
meshgrid(vx,vy)
Computes an (x,y)-grid from the vectors (vx,vy).
For more information, see the MATLAB documentation.

From: https://github.com/ChrisRackauckas/VectorizedRoutines.jl
"""
function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where {T}
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n)
    vy = reshape(vy, m, 1)
    return (repeat(vx, m, 1), repeat(vy, 1, n))
end