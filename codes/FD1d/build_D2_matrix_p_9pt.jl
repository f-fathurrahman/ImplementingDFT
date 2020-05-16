"""
Build second derivative matrix using 9-points centered
finite difference approximation.

# Arguments
- `N::Int64`: number of grid points
- `h::Float64`: spacing between grid points

f_xx = (-9*f[i-4] + 128*f[i-3] - 1008*f[i-2] + 8064*f[i-1] - 14350*f[i+0] + 
        8064*f[i+1]-1008*f[i+2]+128*f[i+3]-9*f[i+4])/(5040*h^2)

Periodic boundary condition
"""
function build_D2_matrix_p_9pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-4
        mat[i,i]   = -14350.0
        mat[i,i+1] =   8064.0
        mat[i,i+2] =  -1008.0
        mat[i,i+3] =    128.0
        mat[i,i+4] =     -9.0
        mat[i+1,i] =   8064.0
        mat[i+2,i] =  -1008.0
        mat[i+3,i] =    128.0
        mat[i+4,i] =     -9.0
    end

    mat[N-3,N-3] = -14350.0
    #
    mat[N-3,N-2] =  8064.0
    mat[N-3,N-1] = -1008.0
    mat[N-3,N]   =   128.0
    #
    mat[N-2,N-3] =  8064.0
    mat[N-1,N-3] = -1008.0
    mat[N,N-3]   =   128.0

    mat[N-2,N-2] = -14350.0
    mat[N-2,N-1] =  8064.0
    mat[N-2,N]   = -1008.0
    mat[N-1,N-2] =  8064.0
    mat[N,N-2]   = -1008.0

    mat[N-1,N-1] = -14350.0
    mat[N-1,N]   = 8064.0
    mat[N,N-1]   = 8064.0

    mat[N,N] = -14350.0

    mat[1,N]   =  8064.0
    mat[1,N-1] = -1008.0
    mat[1,N-2] =   128.0
    mat[1,N-3] =    -9.0
    #
    mat[N  ,1] =  8064.0
    mat[N-1,1] = -1008.0
    mat[N-2,1] =   128.0
    mat[N-3,1] =    -9.0

    mat[2,N  ] = -1008.0
    mat[2,N-1] =   128.0
    mat[2,N-2] =    -9.0
    #
    mat[N  ,2] = -1008.0
    mat[N-1,2] =   128.0
    mat[N-2,2] =    -9.0

    mat[3,N  ] = 128.0
    mat[3,N-1] =  -9.0
    #
    mat[N  ,3] = 128.0
    mat[N-1,3] =  -9.0

    mat[4,N] = -9.0
    mat[N,4] = -9.0

    return mat/(5040*h^2)
end

