function op_lap(H::Hamiltonian, x::Array{Float64,1})
    # we ask for a vector, given that we will consider the vector to be
    # a nx1 matrix
    # TODO: this can be optimized using rfft
    # TODO: we can surther accelerate this using a in-place multiplication

    ytemp = H.kmul.*fft(x)
    return real(ifft(ytemp))
end

function op_lap(H::Hamiltonian, x::Array{Float64,2})
    # application of the laplacian part to a matrix.
    # TODO: this can be optimized using rfft
    # TODO: we can further accelerate this using a in-place multiplication
    ytemp = broadcast(*, H.kmul, fft(x,1))
    return real(ifft(ytemp,1))
end

function op_inv_lap(H::Hamiltonian, x::Array{Float64,1})
    # inverse laplacian, to be used as a preconditioner for the
    # lobpcg algorithm
    inv_kmul = zeros(size(H.kmul))
    inv_kmul[1] = 0
    inv_kmul[2:end] = 1 ./H.kmul[2:end]

    ytemp = inv_kmul.*fft(x)
    return real(ifft(ytemp))
end

function lap_opt(H::Hamiltonian, x::Array{Float64,1})
    # we ask for a 2 vector, given that we will consider the vector to be
    xFourier = rfft(x)
    laplacian_fourier_mult!(xFourier, H.Ls)
    return irfft(xFourier, H.Ns )
end

function lap_opt!(H::Hamiltonian, x::Array{Float64,1})
    # we ask for a 2 vector, given that we will consider the vector to be
    xFourier = rfft(x)
    laplacian_fourier_mult!(xFourier, H.Ls)
    x[:] = irfft(xFourier, H.Ns )
end

function laplacian_fourier_mult!(R::Vector{ComplexF64}, Ls::Float64 )
    c = (2 * pi / Ls)^2/2
    @simd for ii = 1:length(R)
        @inbounds R[ii] = (ii-1)^2*c*R[ii]
    end
end