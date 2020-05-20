"""
f_x = (-1*f[i-1] + 0*f[i+0] + 1*f[i+1])/(2*h**1)
"""
function build_D1_matrix_p_3pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-1
        mat[i,i+1] =  1.0
        mat[i+1,i] = -1.0
    end
    
    mat[1,N] = -1.0
    mat[N,1] =  1.0

    return mat/(2*h)

end