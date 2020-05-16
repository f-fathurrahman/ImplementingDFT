function init_LF1d_c_grid( x_min::Float64, x_max::Float64, N::Int64 )
    L = x_max - x_min
    h = L/(N + 1)
    x = zeros(Float64,N)
    for i = 1:N
        x[i] = x_min + i*L/(N+1)
    end
    return x, h
end

init_LF1d_c_grid( X, N ) = init_LF1d_c_grid( X[1], X[2], N )