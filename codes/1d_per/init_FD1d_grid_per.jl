function init_FD1d_grid_per( x_min::Float64, x_max::Float64, N::Int64 )
    L = x_max - x_min
    h = L/N
    x = zeros(Float64,N)
    for i = 1:N
        x[i] = x_min + (i-1)*h
    end
    return x, h
end

init_FD1d_grid_per( X::Tuple{Float64,Float64}, N::Int64) = init_FD1d_grid_per( X[1], X[2], N )