mutable struct Hamiltonian
    grid
    Laplacian::SparseMatrixCSC{Float64,Int64}
    V_Ps_loc::Vector{Float64}
    V_Hartree::Vector{Float64}
    rhoe::Vector{Float64}
    Nelectrons::Int64
    Focc::Vector{Float64}
    atoms::Atoms
    precKin
    psolver
    gvec::Union{Nothing,GVectors}
end

function Hamiltonian( atoms::Atoms, grid, V_Ps_loc; Nelectrons=2 )

    Laplacian = build_nabla2_matrix( grid )

    # Initialize gvec for periodic case
    if grid.pbc == (true,true,true)
        gvec = GVectors(grid)
    else
        gvec = nothing
    end

    if iseven(Nelectrons)
        Nstates = round(Int64,Nelectrons/2)
        Focc = 2.0*ones(Float64,Nstates)
    else
        Nstates = round(Int64, (Nelectrons+1)/2)
        Focc = 2.0*ones(Float64,Nstates)
        Focc[Nstates] = 1.0
    end

    Npoints = grid.Npoints
    V_Hartree = zeros(Float64, Npoints)
    Rhoe = zeros(Float64, Npoints)

    @printf("Building preconditioners ...")
    precKin = aspreconditioner( ruge_stuben(-0.5*Laplacian) )
    @printf("... done\n")

    if grid.pbc == (false,false,false)
        psolver = PoissonSolverDAGE(grid)
    else
        psolver = PoissonSolverFFT(grid)
    end

    return Hamiltonian( grid, Laplacian, V_Ps_loc, V_Hartree, Rhoe,
        Nelectrons, Focc, atoms, precKin, psolver, gvec )
end

import Base: *
function *( Ham::Hamiltonian, psi::Matrix{Float64} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    Hpsi = zeros(Float64,Nbasis,Nstates)
    Hpsi = -0.5 * Ham.Laplacian * psi
    for ist in 1:Nstates, ip in 1:Nbasis
        Hpsi[ip,ist] = Hpsi[ip,ist] + ( Ham.V_Ps_loc[ip] + Ham.V_Hartree[ip] ) * psi[ip,ist]
    end
    return Hpsi
end

function update!( Ham::Hamiltonian, Rhoe::Vector{Float64} )
    Ham.rhoe[:] = Rhoe[:]
    Ham.V_Hartree = Poisson_solve( Ham.psolver, Ham.grid, Rhoe )
    return
end
