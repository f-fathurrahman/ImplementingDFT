function build_D2_matrix_LF1d_sinc( x::Array{Float64,1}, h::Float64, N::Int64 )
    D2 = zeros(Float64,N,N)
    # Diagonal part
    for i in 1:N
        D2[i,i] = -pi^2 / 3.0 / h^2
    end
    # Off-diagonal
    for j in 1:N
        for i in j+1:N
            D2[i,j] = -2.0*(-1.0)^(i-j)/( x[i] - x[j] )^2
            D2[j,i] = D2[i,j]
        end
    end
    return D2
end