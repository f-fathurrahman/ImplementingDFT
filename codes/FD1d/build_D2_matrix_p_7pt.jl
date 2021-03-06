"""
Build second derivative matrix using 7-points centered
finite difference approximation.

# Arguments
- `N::Int64`: number of grid points
- `h::Float64`: spacing between grid points

f_xx = (2*f[i-3] - 27*f[i-2] + 270*f[i-1] - 490*f[i+0] + 270*f[i+1] - 27*f[i+2] + 2*f[i+3] )/(180*h^2)

Periodic boundary condition
"""
function build_D2_matrix_p_7pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-3
        mat[i,i] = -490.0
        mat[i,i+1] = 270.0
        mat[i,i+2] = -27.0
        mat[i,i+3] = 2.0
        mat[i+1,i] = 270.0
        mat[i+2,i] = -27.0
        mat[i+3,i] = 2.0
    end

    mat[N,N] = -490.0

    mat[N-2,N-2] = -490.0
    #
    mat[N-2,N-1] = 270.0
    mat[N-2,N] = -27.0
    #
    mat[N-1,N-2] = 270.0
    mat[N,N-2] = -27.0

    mat[N-1,N-1] = -490.0
    mat[N-1,N] = 270.0
    mat[N,N-1] = 270.0

    mat[1,N  ] = 270.0
    mat[1,N-1] = -27.0
    mat[1,N-2] = 2.0
    #
    mat[N  ,1] = 270.0
    mat[N-1,1] = -27.0
    mat[N-2,1] = 2.0

    mat[2,N  ] = -27.0
    mat[2,N-1] = 2.0
    #
    mat[N  ,2] = -27.0
    mat[N-1,2] = 2.0

    mat[3,N] = 2.0
    mat[N,3] = 2.0

    return mat/(180.0*h^2)
end

