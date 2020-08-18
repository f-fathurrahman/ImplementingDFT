const AMG_PREC_TYPE = typeof( aspreconditioner(ruge_stuben(speye(1))) )

mutable struct Hamiltonian
    grid::Union{FD3dGrid,LF3dGrid}
    Laplacian::SparseMatrixCSC{Float64,Int64}
    V_Ps_loc::Vector{Float64}
    V_Hartree::Vector{Float64}
    V_XC::Vector{Float64}
    electrons::Electrons
    rhoe::Vector{Float64}
    atoms::Atoms
    precKin::Union{AMG_PREC_TYPE,ILU0Preconditioner}
    psolver::Union{PoissonSolverDAGE,PoissonSolverFFT}
    energies::Energies
    gvec::Union{Nothing,GVectors}
end

"""
Build a Hamiltonian with given FD grid and local potential.
"""
function Hamiltonian( atoms::Atoms, grid, V_Ps_loc;
    Nelectrons=2, Nstates_extra=0,
    stencil_order=9
)
    
    # Need better mechanism for this
    if typeof(grid) == FD3dGrid
        Laplacian = build_nabla2_matrix( grid, stencil_order=stencil_order )
    else
        Laplacian = build_nabla2_matrix( grid )
    end

    # Initialize gvec for periodic case
    if grid.pbc == (true,true,true)
        gvec = GVectors(grid)
    else
        gvec = nothing
    end

    Npoints = grid.Npoints
    V_Hartree = zeros(Float64, Npoints)

    V_XC = zeros(Float64, Npoints)
    Rhoe = zeros(Float64, Npoints)

    @printf("Building preconditioners ...")
    precKin = aspreconditioner( ruge_stuben(-0.5*Laplacian) )
    @printf("... done\n")

    electrons = Electrons( Nelectrons, Nstates_extra=Nstates_extra )

    if grid.pbc == (false,false,false)
        psolver = PoissonSolverDAGE(grid)
    else
        psolver = PoissonSolverFFT(grid)
    end

    energies = Energies()
    
    return Hamiltonian( grid, Laplacian, V_Ps_loc, V_Hartree, V_XC, electrons,
                        Rhoe, atoms, precKin, psolver, energies, gvec )
end

import Base: *
function *( Ham::Hamiltonian, psi::Matrix{Float64} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    Hpsi = zeros(Float64,Nbasis,Nstates)
    Hpsi = -0.5*Ham.Laplacian * psi
    for ist in 1:Nstates, ip in 1:Nbasis
        Hpsi[ip,ist] = Hpsi[ip,ist] + ( Ham.V_Ps_loc[ip] + Ham.V_Hartree[ip] + Ham.V_XC[ip] ) * psi[ip,ist]
    end
    return Hpsi
end

function update!( Ham::Hamiltonian, Rhoe::Vector{Float64} )
    Ham.rhoe = Rhoe
    Ham.V_Hartree = Poisson_solve( Ham.psolver, Ham.grid, Rhoe )
    Ham.V_XC = calc_Vxc_VWN( XCCalculator(), Rhoe )
    return
end
