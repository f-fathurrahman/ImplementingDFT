mutable struct Hamiltonian
    grid::Union{FD3dGrid,LF3dGrid}
    Laplacian::SparseMatrixCSC{Float64,Int64}
    V_Ps_loc::Vector{Float64}
    V_Hartree::Vector{Float64}
    V_XC::Vector{Float64}
    electrons::Electrons
    rhoe::Vector{Float64}
    precKin
    precLaplacian
    energies::Energies
end

"""
Build a Hamiltonian with given FD grid and local potential.
"""
function Hamiltonian( grid, ps_loc_func::Function;
    Nelectrons=2, Nstates_extra=0,
    func_1d=build_D2_matrix_5pt
)
    
    # Need better mechanism for this
    if typeof(grid) == FD3dGrid
        Laplacian = build_nabla2_matrix( grid, func_1d=func_1d )
    else
        Laplacian = build_nabla2_matrix( grid )
    end
    V_Ps_loc = ps_loc_func( grid )

    Npoints = grid.Npoints
    V_Hartree = zeros(Float64, Npoints)

    V_XC = zeros(Float64, Npoints)
    Rhoe = zeros(Float64, Npoints)

    @printf("Building preconditioners ...")
    precKin = aspreconditioner( ruge_stuben(-0.5*Laplacian) )
    precLaplacian = aspreconditioner( ruge_stuben(Laplacian) )
    @printf("... done\n")

    electrons = Electrons( Nelectrons, Nstates_extra=Nstates_extra )

    energies = Energies()
    return Hamiltonian( grid, Laplacian, V_Ps_loc, V_Hartree, V_XC, electrons,
                        Rhoe, precKin, precLaplacian, energies )
end


function Hamiltonian( grid, V_loc_func::Array{Float64,1};
    Nelectrons=2, Nstates_extra=0,
    func_1d=build_D2_matrix_5pt
)

    Laplacian = build_nabla2_matrix( grid, func_1d=func_1d )
    V_Ps_loc = V_loc_func

    Npoints = grid.Npoints
    V_Hartree = zeros(Float64, Npoints)

    V_XC = zeros(Float64, Npoints)
    Rhoe = zeros(Float64, Npoints)

    @printf("Building preconditioners ...")
    precKin = aspreconditioner( ruge_stuben(-0.5*Laplacian) )
    precLaplacian = aspreconditioner( ruge_stuben(Laplacian) )
    @printf("... done\n")

    electrons = Electrons( Nelectrons, Nstates_extra=Nstates_extra )

    energies = Energies()
    return Hamiltonian( grid, Laplacian, V_Ps_loc, V_Hartree, V_XC, electrons,
                        Rhoe, precKin, precLaplacian, energies )
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
    Ham.V_Hartree = Poisson_solve_PCG( Ham.Laplacian, Ham.precLaplacian, Rhoe, 1000, verbose=false, TOL=1e-10 )
    Ham.V_XC = excVWN( Rhoe ) + Rhoe .* excpVWN( Rhoe )
    return
end
