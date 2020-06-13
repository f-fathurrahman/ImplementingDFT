function init_LF1d_p_grid( x_min::Float64, x_max::Float64, N::Int64 )
    
    @assert !iseven(N)

    L = x_max - x_min
    x = zeros(Float64,N)
    for i in 1:N
        x[i] = x_min + 0.5*L*(2*i-1)/N
    end
    h = x[2] - x[1]
    return x, h
end

init_LF1d_p_grid( X, N ) = init_LF1d_p_grid( X[1], X[2], N )