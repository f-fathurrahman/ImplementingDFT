function init_LF1d_sinc_grid(h::Float64, N::Int64)
    # Choice for x_min and x_max
    x_min = -(N-1)*h/2.0
    #
    x = zeros(Float64,N)
    for i in 1:N
        x[i] = x_min + (i-1)*h
    end
    return x
end

function init_LF1d_sinc_grid( x_min::Float64, x_max::Float64, N::Int64 )
    @assert x_min < 0.0
    @assert x_max > 0.0
    @assert abs(x_min) ≈ x_max
    h = 2*x_max/(N-1)
    x = zeros(Float64,N)
    for i in 1:N
        x[i] = x_min + (i-1)*h
    end
    return x, h
end

init_LF1d_sinc_grid( X, N ) = init_LF1d_sinc_grid( X[1], X[2], N )