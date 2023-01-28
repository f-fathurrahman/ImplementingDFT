mutable struct Potentials
    Ions::Array{Float64,1}
    Hartree::Array{Float64,1}
    XC::Array{Float64,2}  # spin dependent
    Total::Array{Float64,2}
end

mutable struct Hamiltonian1d
    atoms::Atoms1d
    grid::FD1dGrid
    electrons::Electrons
    rhoe::Matrix{Float64}
    potentials::Potentials
    Kmat::Matrix{Float64}
    xc_calc::LibxcXCCalculator
    energies::Energies
end


# XXX Make this type parametric according to the basis used:
# FD1dGrid, LF1dGrid, or PW1dGrid
function Hamiltonian1d(
    atoms::Atoms1d,
    Npoints::Int64;
    basis=:fd,
    Nstates_extra=0,
    Nspin=1,
)

    # Limitations on the current implementation
    @assert Nspin == 1
    @assert atoms.pbc == false

    #if basis == :fd
    xmin = -atoms.L/2
    xmax =  atoms.L/2
    grid = FD1dGrid( (xmin, xmax), Npoints, pbc=atoms.pbc)
    Kmat = -0.5*build_D2_matrix_9pt(Npoints, grid.hx)

    #elseif basis_type == :pw
    #    gvec = GVectors1d(L, Ns)
    #end

    hx = grid.hx

    electrons = Electrons(atoms, Nstates_extra=Nstates_extra, Nspin=Nspin)

    rhoe = ones(Float64, Npoints, Nspin)
    @views rhoe[:] = rhoe[:] / ( sum(rhoe)*hx ) * electrons.Nelectrons
    println("integ rhoe = ", sum(rhoe)*hx)
    
    potentials = Potentials(
        zeros(Float64, Npoints),
        zeros(Float64, Npoints),
        zeros(Float64, Npoints, Nspin),
        zeros(Float64, Npoints, Nspin)
    )

    xc_calc = LibxcXCCalculator()

    energies = Energies()

    return Hamiltonian1d( atoms, grid, electrons, rhoe, potentials, Kmat, xc_calc, energies )
end

# This function can be optimized
function op_H(Ham::Hamiltonian1d, psi::Union{Matrix{Float64},Vector{Float64}})
    Kmat = Ham.Kmat
    Vtot = Ham.potentials.Total
    Hmat = Kmat + diagm( 0 => Vtot[:,1] )
    return Hmat * psi
end

import Base.*
*(Ham::Hamiltonian1d, psi) = op_H(Ham, psi)


function prec_invK!(Ham::Hamiltonian1d, v, Kv)
    Kv[:,:] = inv(Ham.Kmat)*v
    return
end

# Ideal preconditioner (expensive to calculate)
function prec_invHam!(Ham::Hamiltonian1d, v, Kv)
    Kmat = Ham.Kmat
    Vtot = Ham.potentials.Total
    Hmat = Kmat + diagm( 0 => Vtot[:,1] )
    λ = eigvals(Hmat)
    fill!(Kv, 0.0)
    for i in 1:size(v,2)
        @views Kv[:,i] = inv(Hmat - λ[i]*I)*v[:,i]
    end
    return
end


# For plane wave basis
#function get_matrix( Ham::Hamiltonian1d )
#    Gkin = 0.5*Ham.gvec.G2
#    Ng = Ham.gvec.Ng
#    # create the matrix version of the Hamiltonian
#    Hmat = real( ifft( diagm(0 => Gkin) * fft(Matrix{Float64}(I, Ng, Ng ), 1), 1) )
#    Hmat += diagm( 0 => Ham.Vtotal[:] )
#    # symmetrize
#    return 0.5*(Hmat + Hmat')
#end


import Base: show
function show( io::IO, Ham::Hamiltonian1d; header=true )
    if header
        @printf("\n")
        @printf("                                -------------\n")
        @printf("                                Hamiltonian1d\n")
        @printf("                                -------------\n")
        @printf("\n")
    end
    @printf(io, "size (MiB) = %18.5f\n", Base.summarysize(Ham)/1024/1024)
    println(io, "")
    show(io, Ham.atoms)
    show(io, Ham.electrons)
end
show( Ham::Hamiltonian1d; header=true ) = show( stdout, Ham, header=header )