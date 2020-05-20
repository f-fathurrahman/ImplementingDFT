"""
f_x = (-1*f[i-3] + 9*f[i-2] - 45*f[i-1] + 45*f[i+1] - 9*f[i+2] + 1*f[i+3])/(60*h**1)
"""
function build_D1_matrix_7pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-3
        #
        mat[i,i+1] = 45.0
        mat[i,i+2] = -9.0
        mat[i,i+3] =  1.0
        #
        mat[i+1,i] = -45.0
        mat[i+2,i] =   9.0
        mat[i+3,i] =  -1.0
    end

    mat[N  ,N-1] = -45.0
    mat[N  ,N-2] =   9.0
    #
    mat[N-2,N  ] =  -9.0
    mat[N-1,N  ] =  45.0

    mat[N-1,N-2] = -45.0
    #
    mat[N-2,N-1] =  45.0

    return mat/(60*h)

end