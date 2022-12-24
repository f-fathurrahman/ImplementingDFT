using Printf
using LaTeXStrings

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("INC_sch_1d.jl")

function my_gaussian(x; α=1.0)
    return exp(-α*x^2)
end


x0 = 1.0 # center of Gaussian
L = 5.0 # periodicity or lattice size
A = 0.0
B = A + L # two times of Lattice
N = 50

xgrid, dx = init_FD1d_p_grid( A, B, N )

fx = zeros(Float64,N)
# Use nearest neighbors
for nn in [-1,0,1]
    fx += my_gaussian.(xgrid .- x0 .+ nn*L, α=2.0)
end

plt.clf()
plt.plot(xgrid, fx, label="using nearest cell")

fx2 = zeros(Float64, N)
fx2[:] = fx[:] # copy
# Use additional neighbors
for nn in [-2, 2]
    fx2 += my_gaussian.(xgrid .- x0 .+ nn*L, α=2.0)
end

plt.plot(xgrid, fx2, label="using 2nd-nearest cell")

plt.legend()
plt.grid(true)
plt.savefig("IMG_periodic_gauss.png", dpi=150)


plt.clf()
plt.plot(xgrid, fx - fx2, label="fx-fx2")
plt.legend()
plt.grid(true)
plt.savefig("IMG_diff_fx.png", dpi=150)

