"""
Build second derivative matrix using 3-points centered
finite difference approximation.

# Arguments
- `N::Int64`: number of grid points
- `h::Float64`: spacing between grid points

Periodic version
"""
function build_D2_matrix_3pt_per( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-1
        mat[i,i] = -2.0
        mat[i,i+1] = 1.0
        mat[i+1,i] = mat[i,i+1]
    end
    mat[1,N] =  1.0
    mat[N,N] = -2.0
    mat[N,1] =  1.0
    return mat/h^2
end

