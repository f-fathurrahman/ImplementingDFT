"""
f_x = (-2*f[i-5] + 25*f[i-4] - 150*f[i-3] + 600*f[i-2] - 2100*f[i-1] + 
       2100*f[i+1] - 600*f[i+2] + 150*f[i+3] - 25*f[i+4] + 2*f[i+5])/(2520*h**1)
"""
function build_D1_matrix_p_11pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-5
        #
        mat[i,i+1] = 2100.0
        mat[i,i+2] = -600.0
        mat[i,i+3] =  150.0
        mat[i,i+4] =  -25.0
        mat[i,i+5] =    2.0
        #
        mat[i+1,i] = -2100.0
        mat[i+2,i] =   600.0
        mat[i+3,i] =  -150.0
        mat[i+4,i] =    25.0
        mat[i+5,i] =    -2.0
    end

    mat[N  ,N-1] = -2100.0
    mat[N  ,N-2] =   600.0
    mat[N  ,N-3] =  -150.0
    mat[N  ,N-4] =    25.0
    #
    mat[N-1,N  ] =  2100.0
    mat[N-2,N  ] =  -600.0
    mat[N-3,N  ] =   150.0
    mat[N-4,N  ] =   -25.0


    mat[N-1,N-2] = -2100.0
    mat[N-1,N-3] =   600.0
    mat[N-1,N-4] =  -150.0
    #
    mat[N-2,N-1] =  2100.0
    mat[N-3,N-1] =  -600.0
    mat[N-4,N-1] =   150.0


    mat[N-2,N-3] = -2100.0
    mat[N-2,N-4] =   600.0
    #
    mat[N-3,N-2] =  2100.0
    mat[N-4,N-2] =  -600.0


    mat[N-3,N-4] = -2100.0
    mat[N-4,N-3] =  2100.0

    # Periodic

    mat[1  ,N  ] = -2100.0
    mat[1  ,N-1] =   600.0
    mat[1  ,N-2] =  -150.0
    mat[1  ,N-3] =    25.0
    mat[1  ,N-4] =    -2.0
    #
    mat[N  ,  1] =  2100.0
    mat[N-1,  1] =  -600.0
    mat[N-2,  1] =   150.0
    mat[N-3,  1] =   -25.0
    mat[N-4,  1] =     2.0

    mat[2  ,N  ] =   600.0
    mat[2  ,N-1] =  -150.0
    mat[2  ,N-2] =    25.0
    mat[2  ,N-3] =    -2.0
    #
    mat[N  ,  2] =  -600.0
    mat[N-1,  2] =   150.0
    mat[N-2,  2] =   -25.0
    mat[N-3,  2] =     2.0


    mat[3  ,N  ] =  -150.0
    mat[3  ,N-1] =    25.0
    mat[3  ,N-2] =    -2.0
    #
    mat[N  ,  3] =   150.0
    mat[N-1,  3] =   -25.0
    mat[N-2,  3] =     2.0

    mat[4  ,N  ] =    25.0
    mat[4  ,N-1] =    -2.0
    #
    mat[N  ,  4] =    -25.0
    mat[N-1,  4] =     2.0

    mat[5,N] = -2.0
    mat[N,5] =  2.0

    return mat/(2520*h)

end