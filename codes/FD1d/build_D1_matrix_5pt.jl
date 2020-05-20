"""
f_x = (1*f[i-2] - 8*f[i-1] + 8*f[i+1] - 1*f[i+2])/(12*h**1)
"""
function build_D1_matrix_5pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-2
        #
        mat[i,i+1] =  8.0
        mat[i,i+2] = -1.0
        #
        mat[i+1,i] = -8.0
        mat[i+2,i] =  1.0
    end
    
    mat[N,  N-1] = -8.0
    mat[N-1,N  ] =  8.0

    return mat/(12*h)

end