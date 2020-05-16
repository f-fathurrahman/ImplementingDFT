function build_D2_matrix_LF1d_c( L::Float64, N::Int64 )
    
    D2 = zeros(Float64,N,N)
    
    # Diagonal part
    pre = -pi^2/(2.0*L^2)
    for i = 1:N
        t1 = ( 2.0*(N+1)^2 + 1 )/3.0
        t2 = sin( i*pi/(N+1) )^2
        D2[i,i] = pre*( t1 - 1.0/t2 )
    end
    
    # Off-diagonal
    for l = 1 : N
        for j = l+1 : N
            nnm = l - j
            nnp = l + j
            pre = -pi^2 / (2*L^2) * (-1.0)^nnm
            t1 = sin( pi*nnm/2.0/(N+1) )^2
            t2 = sin( pi*nnp/2.0/(N+1) )^2
            #
            D2[l,j] = pre*( 1.0/t1 - 1.0/t2 )
            D2[j,l] = D2[l,j]
        end
    end
    
    return D2
end

build_D2_matrix_LF1d_c(x_min, x_max, N) = build_D2_matrix_LF1d_c( x_max - x_min, N )
