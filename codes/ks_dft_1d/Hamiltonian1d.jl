mutable struct Hamiltonian1d
    atoms::Atoms1d
    electrons::Electrons
    rhoe::Matrix{Float64}
    rhoa::Vector{Float64}
    drhoa::Matrix{Float64}
    VHartree::Vector{Float64}
    Vtotal::Vector{Float64}
end

function Hamiltonian1d(atoms, dx_in, κ, ε0; Nstates_extra=0, Nspin=1)

    @assert Nspin == 1 # current implementation

    # Setup grid
    L = atoms.L
    Ns = round(Int64, L/dx_in)
    dx = L/Ns
    Npoints = Ns

    gvec = GVectors1d(L, Ns)
    electrons = Electrons(atoms, Nstates_extra=Nstates_extra, Nspin=Nspin)

    rhoa, drhoa = init_pseudocharges(atoms, dx, Ns)

    println("integ rhoa = ", sum(rhoa)*dx)

    rhoe = zeros(Float64, Ns, Nspin)
    @views rhoe[:,1] = -rhoa
    rhoe[:] = rhoe[:] / (sum(rhoe)*dx) * electrons.Nelectrons
    println("integ rhoe = ", sum(rhoe)*dx)
    
    VHartree = zeros(Float64, Npoints)
    Vtotal = zeros(Float64, Npoints)

    return Hamiltonian1d(atoms, gvec, Ns, κ, ε0, dx, electrons,
        rhoe, rhoa, drhoa, VHartree, Vtotal)
end


#
# Create the pseudo-charge for Coulomb interaction
# $\rho_a = \sum_{N_{\text{atoms}}} - Z_j/\sqrt{2 \pi \sigma_i}$
#
function init_pseudocharges( atoms::Atoms1d, dx::Float64, Ns::Int64 )

    Natoms = atoms.Natoms
    Zvals = atoms.Zvals
    σ = atoms.σ
    atpos = atoms.positions
    L = atoms.L
    Npoints = Ns

    # rgrid[ip] = (ip-1)*dx

    # Calculate the total pseudo-charge
    # NOTE: The pseudo-charge should not extend over twice the domain size!
    rhoa = zeros(Float64, Npoints)
    for ia in 1:Natoms
        for ip in 1:Npoints
            d = atpos[ia] - (ip-1)*dx
            d = d - round(Int64, d/L)*L # why? (PBC)
            # Note: Z has minus sign in front of it!
            rhoa[ip] += -Zvals[ia]/sqrt(2π*σ[ia]^2) * exp(-0.5*(d/σ[ia])^2)
        end
    end

    drhoa = zeros(Float64, Npoints, Natoms) # FIXME CHECK THIS
    for ia in 1:Natoms
        for ip in 1:Natoms
            d = atpos[ia] - (ip-1)*dx
            d = d - round(Int64, d/L)*L
            drhoa[ip,ia] += -Zvals[ia]/sqrt(2π*σ[ia]^2)/σ[ia]^2 * exp(-0.5*(d/σ[ia])^2)*d
        end
    end
    return rhoa, drhoa
end


function get_matrix( Ham::Hamiltonian1d )
    Gkin = 0.5*Ham.gvec.G2
    Ng = Ham.gvec.Ng
    # create the matrix version of the Hmailtonian
    Hmat = real( ifft( diagm(0 => Gkin) * fft(Matrix{Float64}(I, Ng, Ng ), 1), 1) )
    Hmat += diagm( 0 => Ham.Vtotal[:] )
    # symmetrize
    return 0.5*(Hmat + Hmat')
end


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