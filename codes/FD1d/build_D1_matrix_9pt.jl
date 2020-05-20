"""
f_x = (3*f[i-4] - 32*f[i-3] + 168*f[i-2] - 672*f[i-1] + 
       672*f[i+1] - 168*f[i+2] + 32*f[i+3]- 3*f[i+4] ) / ( 840*h**1 )
"""
function build_D1_matrix_9pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-4
        #
        mat[i,i+1] =  672.0
        mat[i,i+2] = -168.0
        mat[i,i+3] =   32.0
        mat[i,i+4] =   -3.0
        #
        mat[i+1,i] = -672.0
        mat[i+2,i] =  168.0
        mat[i+3,i] =  -32.0
        mat[i+4,i] =    3.0
    end

    mat[N  ,N-1] = -672.0
    mat[N  ,N-2] =  168.0
    mat[N  ,N-3] =  -32.0
    #
    mat[N-1,N  ] =  672.0
    mat[N-2,N  ] = -168.0
    mat[N-3,N  ] =   32.0

    mat[N-1,N-2] = -672.0
    mat[N-1,N-3] =  168.0
    #
    mat[N-2,N-1] =  672.0
    mat[N-3,N-1] = -168.0

    mat[N-2,N-3] = -672.0
    #
    mat[N-3,N-2] =  672.0

    return mat/(840*h)

end