"""
Build second derivative matrix using 5-points centered
finite difference approximation.

# Arguments
- `N::Int64`: number of grid points
- `h::Float64`: spacing between grid points

f_xx = (-1*f[i-2] + 16*f[i-1]- 30*f[i+0] + 16*f[i+1] - 1*f[i+2] )/(12*h^2)
"""
function build_D2_matrix_5pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-2
        mat[i,i] = -30.0
        mat[i,i+1] = 16.0
        mat[i,i+2] = -1.0
        mat[i+1,i] = 16.0
        mat[i+2,i] = -1.0
    end

    mat[N-1,N-1] = -30.0
    mat[N-1,N] = 16.0
    mat[N,N-1] = 16.0

    mat[N,N] = -30.0

    return mat/(12*h^2)
end