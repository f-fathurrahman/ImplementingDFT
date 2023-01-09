include("test_direct_min_01.jl")

Ham = init_Hamiltonian()

hx = Ham.grid.hx
Npoints = Ham.grid.Npoints
Nstates = Ham.electrons.Nstates

psi = ortho_sqrt( rand(Float64, Npoints, Nstates) )
psi[:] = psi[:]/sqrt(hx)
E0 = calc_KohnSham_Etotal!(Ham, psi)
g0 = calc_grad(Ham, psi)

Δ = 1e-10
dW = randn(Float64, Npoints, Nstates)
psi1 = ortho_sqrt( psi + Δ*dW )
psi1[:] = psi1[:]/sqrt(hx)
E1 = calc_KohnSham_Etotal!(Ham, psi1)

dE = 2*real( sum( g0 .* (Δ*dW) ) )*hx

println("dE = ", dE)
println("E1 - E0 = ", E1 - E0)
println("ratio (should be close to 1) = ", abs((E1 - E0)/dE))