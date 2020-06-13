function build_D2_matrix_LF1d_p( L::Float64, N::Int64 )
    
    @assert !iseven(N)

    D2 = zeros(Float64,N,N)    
    Nprimed = round(Int64,(N-1)/2)

    # Diagonal elements
    for j in  1:N
        D2[j,j] = -( 2*pi/L )^2 * Nprimed * (Nprimed + 1)/3.0
    end
    
    # Off diagonal elements
    for j in 1:N, l = (j+1):N
        n = j - l        
        tt3 = (2.0*pi/L)^2 * (-1.0)^n * cos(pi*n/N)
        tt4 = 2.0*sin(pi*n/N)^2        
        D2[j,l] = -tt3/tt4
        D2[l,j] = -tt3/tt4
    end

    return D2
end

build_D2_matrix_LF1d_p(x_min, x_max, N) = build_D2_matrix_LF1d_p( x_max - x_min, N )
