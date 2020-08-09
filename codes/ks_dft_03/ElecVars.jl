#
# Subspace rotation matrices
#

mutable struct SubspaceRotations
    prev::Matrix{Float64}
    prevC::Matrix{Float64}
    prevCinv::Matrix{Float64}
end

function SubspaceRotations( Nstates::Int64 )
    prev[i] = diagm( 0 => ones(Float64,Nstates) )
    prevC[i] = diagm( 0 => ones(Float64,Nstates) )
    prevCinv[i] = diagm( 0 => ones(Float64,Nstates) )
    return SubspaceRotations(prev, prevC, prevCinv)
end


#
# Electronic gradients
#

mutable struct ElecGradient
    psi::Matrix{Float64}
    Haux::Matrix{Float64}
end

function ElecGradient( Ham )
    Nstates = Ham.electrons.Nstates
    Nbasis = Ham.grid.Npoints
    psi = zeros(Float64, Nbasis, Nstates)
    Haux = zeros(ComplexF64, Nstates, Nstates)
    return ElecGradient(psi, Haux)
end

import Base: -
function -(e1::ElecGradient, e2::ElecGradient)
    return ElecGradient( e1.psi - e2.psi, e1.Haux - e2.Haux )
end

#
# Electronic variables (wave functions and subspace Hamiltonian)
#

mutable struct ElecVars
    psi::Matrix{Float64}
    Hsub::Matrix{Float64}
    Hsub_eigs::Vector{Float64}
end

function ElecVars( Ham::Hamiltonian )
    psi = rand( Float64, Ham.grid.Npoints, Ham.electrons.Nstates )
    ortho_sqrt!( psi, Ham.grid.dVol )
    return ElecVars( Ham, psi )
end

function ElecVars( Ham::Hamiltonian, psi::Matrix{Float64} )

    Rhoe = zeros(Float64, Ham.grid.Npoints)
    calc_rhoe!(Ham, psi, Rhoe)
    update!(Ham, Rhoe)    
    
    Hsub = psi' * op_H(Ham, psi) * Ham.grid.dVol
    Hsub_eigs = eigvals(Hermitian(Hsub))
    return ElecVars( psi, Hsub, Hsub_eigs )
end

import Base: show
function show( io::IO, evars::ElecVars )
    #
    @printf(io, "Hsub = \n")
    Nstates = size(evars.Hsub,1)
    for i in 1:Nstates
        for j in 1:Nstates
            @printf(io, "%10.5f ", evars.Hsub[i,j])
        end
        @printf(io, "\n")
    end
    #
    @printf(io, "Hsub_eigs = \n")
    for i in 1:Nstates
        @printf(io, "%4d %10.5f\n", i, evars.Hsub_eigs[i])
    end
end
show( evars::ElecVars ) = show(stdout, evars)
