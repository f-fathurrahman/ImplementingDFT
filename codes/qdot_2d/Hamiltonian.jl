const AMG_PREC_TYPE = typeof( aspreconditioner(ruge_stuben(speye(1))) )

mutable struct Hamiltonian
    grid::Union{FD2dGrid,LF2dGrid}
    Laplacian::SparseMatrixCSC{Float64,Int64}
    V_Ps_loc::Vector{Float64}
    V_Hartree::Vector{Float64}
    V_XC::Vector{Float64}
    electrons::Electrons
    rhoe::Vector{Float64}
    precKin::Union{AMG_PREC_TYPE,ILU0Preconditioner}
    energies::Energies
end

function Hamiltonian( grid, V_loc::Array{Float64,1};
    Nelectrons=2, Nstates_extra=0,
    stencil_order=9, prec_type=:ILU0, verbose=false
)

    # Need better mechanism for this
    if typeof(grid) == FD2dGrid
        Laplacian = build_nabla2_matrix( grid, stencil_order=stencil_order )
    else
        Laplacian = build_nabla2_matrix( grid )
    end

    V_Ps_loc = V_loc

    Npoints = grid.Npoints
    V_Hartree = zeros(Float64, Npoints)

    V_XC = zeros(Float64, Npoints)
    Rhoe = zeros(Float64, Npoints)

    verbose && @printf("Building preconditioners ...")
    if prec_type == :amg
        precKin = aspreconditioner(ruge_stuben(-0.5*Laplacian))
    else
        precKin = ILU0Preconditioner(-0.5*Laplacian)
    end
    verbose && @printf("... done\n")

    electrons = Electrons( Nelectrons, Nstates_extra=Nstates_extra )

    energies = Energies()
    return Hamiltonian( grid, Laplacian, V_Ps_loc, V_Hartree, V_XC, electrons,
                        Rhoe, precKin, energies )
end

# In-place version
function op_H!( Ham::Hamiltonian, psi, Hpsi )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    Hpsi[:,:] = -0.5*Ham.Laplacian * psi
    for ist in 1:Nstates, ip in 1:Nbasis
        Hpsi[ip,ist] = Hpsi[ip,ist] + ( Ham.V_Ps_loc[ip] + Ham.V_Hartree[ip] + Ham.V_XC[ip] ) * psi[ip,ist]
    end
    return
end

function op_H( Ham::Hamiltonian, psi )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    Hpsi = -0.5*Ham.Laplacian * psi
    for ist in 1:Nstates, ip in 1:Nbasis
        Hpsi[ip,ist] = Hpsi[ip,ist] + ( Ham.V_Ps_loc[ip] + Ham.V_Hartree[ip] + Ham.V_XC[ip] ) * psi[ip,ist]
    end
    return Hpsi
end

import Base: *
function *(Ham, psi)
    return op_H(Ham, psi)
end

function update!( Ham::Hamiltonian, Rhoe::Vector{Float64} )
    Ham.rhoe = Rhoe
    Ham.V_Hartree = Poisson_solve_sum( Ham.grid, Rhoe )
    Ham.V_XC = calc_V_xc_2d(Rhoe)
    return
end
