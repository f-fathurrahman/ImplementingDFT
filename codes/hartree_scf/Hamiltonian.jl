mutable struct Hamiltonian
    grid
    Laplacian::SparseMatrixCSC{Float64,Int64}
    V_Ps_loc::Vector{Float64}
    V_Hartree::Vector{Float64}
    rhoe::Vector{Float64}
    precKin
    precLaplacian
end

function Hamiltonian( grid, ps_loc_func::Function )
    Laplacian = build_nabla2_matrix( grid )
    V_Ps_loc = ps_loc_func( grid )
    Npoints = grid.Npoints
    V_Hartree = zeros(Float64, Npoints)

    Rhoe = zeros(Float64, Npoints)

    @printf("Building preconditioners ...")
    precKin = aspreconditioner( ruge_stuben(-0.5*Laplacian) )
    precLaplacian = aspreconditioner( ruge_stuben(Laplacian) )
    @printf("... done\n")
    return Hamiltonian( grid, Laplacian, V_Ps_loc, V_Hartree, Rhoe, precKin, precLaplacian )
end

import Base: *
function *( Ham::Hamiltonian, psi::Matrix{Float64} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    Hpsi = zeros(Float64,Nbasis,Nstates)
    # include occupation number factor
    Hpsi = -0.5*Ham.Laplacian * psi #* 2.0
    for ist in 1:Nstates, ip in 1:Nbasis
        Hpsi[ip,ist] = Hpsi[ip,ist] + ( Ham.V_Ps_loc[ip] + Ham.V_Hartree[ip] ) * psi[ip,ist]
    end
    return Hpsi
end

function update!( Ham::Hamiltonian, Rhoe::Vector{Float64} )
    Ham.rhoe[:] = Rhoe[:]
    Ham.V_Hartree = Poisson_solve_PCG( Ham.Laplacian, Ham.precLaplacian, Rhoe, 1000, verbose=false, TOL=1e-10 )
    return
end
