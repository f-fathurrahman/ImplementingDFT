function build_D1_matrix_LF1d_p( L::Float64, N::Int64 )
    
    @assert !iseven(N)

    D1 = zeros(Float64,N,N)    
    Nprimed = round(Int64,(N-1)/2)

    # Diagonal elements are already zero

    # Off diagonal elements
    for j in 1:N, l in (j+1):N
        n = j - l        
        tt1 = pi/L * (-1.0)^n
        tt2 = sin(pi*n/N)
        D1[j,l] =  tt1/tt2
        D1[l,j] = -tt1/tt2
    end

    return D1
end

build_D1_matrix_LF1d_p(x_min, x_max, N) = build_D1_matrix_LF1d_p( x_max - x_min, N )
