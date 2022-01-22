using Printf

include("INC_sch_1d.jl")

function my_gaussian(x; α=1.0)
    return exp(-α*x^2)
end

function do_integrate(N::Int64)
    A = -5.0
    B = 5.0
    α = 1.0
    xgrid, dx = init_FD1d_grid( A, B, N )
    fx = my_gaussian.(xgrid, α=α)
    exact_res = sqrt(π)/sqrt(α)
    num_res = sum(fx)*dx
    @printf("%5d %18.10f %10.5e\n", N, num_res, abs(exact_res - num_res))
end

for N in [10, 20, 30, 40, 50]
    do_integrate(N)
end
