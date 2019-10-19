mutable struct Hamiltonian
    fdgrid::FD3dGrid
    Laplacian::SparseMatrixCSC{Float64,Int64}
    V_Ps_loc::Vector{Float64}
    V_Hartree::Vector{Float64}
    rhoe::Vector{Float64}
    precKin
    precLaplacian
end

function Hamiltonian( fdgrid::FD3dGrid, ps_loc_func::Function; func_1d=build_D2_matrix_5pt )
    Laplacian = build_nabla2_matrix( fdgrid, func_1d=func_1d )
    V_Ps_loc = ps_loc_func( fdgrid )
    Npoints = fdgrid.Npoints
    V_Hartree = zeros(Float64, Npoints)

    Rhoe = zeros(Float64, Npoints)

    @printf("Building preconditioners ...")
    precKin = aspreconditioner( ruge_stuben(-0.5*Laplacian) )
    precLaplacian = aspreconditioner( ruge_stuben(Laplacian) )
    @printf("... done\n")
    return Hamiltonian( fdgrid, Laplacian, V_Ps_loc, V_Hartree, Rhoe, precKin, precLaplacian )
end

import Base: *
function *( Ham::Hamiltonian, psi::Matrix{Float64} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    Hpsi = zeros(Float64,Nbasis,Nstates)
    # include occupation number factor
    Hpsi = -0.5*Ham.Laplacian * psi * 2.0
    for ist in 1:Nstates, ip in 1:Nbasis
        Hpsi[ip,ist] = Hpsi[ip,ist] + 2.0*( Ham.V_Ps_loc[ip] + Ham.V_Hartree[ip] ) * psi[ip,ist]
    end
    return Hpsi
end

function update!( Ham::Hamiltonian, Rhoe::Vector{Float64} )
    Ham.rhoe[:] = Rhoe[:]
    Ham.V_Hartree = Poisson_solve_PCG( Ham.Laplacian, Ham.precLaplacian, -4*pi*Rhoe, 1000, verbose=false, TOL=1e-10 )
    V_loc = Ham.V_Hartree + Ham.V_Ps_loc
    println("sum(V_Hartree)*dVol = ", sum(Ham.V_Hartree)*Ham.fdgrid.dVol)
    println("sum(V_loc)*dVol = ", sum(V_loc)*Ham.fdgrid.dVol)
    return
end
