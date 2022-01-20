using Printf

include("lobpcg_sep.jl")

mutable struct Hamiltonian
    Ns::Int64
    Ls::Float64
    kmul::Array{Float64,2} # wait until everything is clearer
    dx::Float64
    gridpos::Array{Float64,2}
    posstart::Int64
    posidx::Int64
    rhoa::Array{Float64,2}               # pseudo-charge for atoms (negative)
    rho::Array{Float64,2}                # electron density
    Vhar::Array{Float64,2}               # Hartree potential for both electron and nuclei
    Vtot::Array{Float64,2}               # total potential
    drhoa::Array{Float64,2}              # derivative of the pseudo-charge (on each atom)
    ev::Array{Float64,1}
    psi::Array{Float64,2}
    fermi::Float64
    occ
    nspin::Int64
    Neigs::Int64    # QUESTION: Number of eigenvalues?
    atoms::Atoms
    Eband::Float64              # Total band energy
    Fband::Float64               # Helmholtz band energy
    Ftot::Float64                # Total Helmholtz energy
    YukawaK::Float64             # shift for the potential
    epsil0::Float64
    Tbeta::Float64               # temperature 1beta = 1/T
end

function Hamiltonian(Lat, Nunit, n_extra, dx, atoms,YukawaK, epsil0,Tbeta)
    Ls = Nunit * Lat
    Ns = round(Int64, Ls / dx)
    #
    dx = Ls / Ns
    # defining the grid
    gridpos = zeros(Ns,1) # allocating as a 2D Array
    gridpos[:,1] = transpose(collect(0:Ns-1)).*dx
    posstart = 0
    posidx   = 0
    # Initialize the atom positions
    atoms   = atoms
    Neigs = sum(atoms.nocc) + n_extra

    println("sum(atoms.nocc) = ", sum(atoms.nocc))
    println("Neigs = ", Neigs)

    Ls_glb = Ls
    Ns_glb = Ns

    # we define the Fourier multipliers as an 2D array
    kx = zeros(Ns,1)
    kx[:,1] = vcat( collect(0:Ns/2-1), collect( -Ns/2:-1) )* 2 * pi / Ls
    kmul = kx.^2/2

    rhoa, drhoa = pseudocharge(gridpos, Ls_glb, atoms, YukawaK, epsil0)

    # TODO: we need to figure out the type of each of the fields to properlu
    # initialize them
    rho = zeros(1,1)
    Vhar = zeros(1,1)
    Vtot = zeros(1,1)
    # drhoa = []
    ev = []
    psi = zeros(1,1)
    fermi = 0.0
    occ = []
    Eband = 0.0
    Fband = 0.0
    Ftot = 0.0
    nspin = 1

    return Hamiltonian(Ns, Ls, kmul, dx, gridpos, posstart, posidx, rhoa,
            rho, Vhar, Vtot, drhoa, ev, psi, fermi, occ,nspin, Neigs, atoms,
            Eband, Fband, Ftot, YukawaK, epsil0, Tbeta)
end


# importing the necessary functions to comply with Julia duck typing
using LinearAlgebra
import Base.*
import LinearAlgebra.mul!
import Base.eltype
import Base.size
import LinearAlgebra.issymmetric

function *(H::Hamiltonian, x::Array{Float64,1})
    y_lap  = op_lap(H, x)
    y_vtot = op_Vtot(H, x)
    return y_lap + y_vtot
end

# TODO: this should be optimized in order to take advantage of BLAS 3 operations
function *(H::Hamiltonian, X::Array{Float64,2})
    # Function  to  overload the application of the Hamiltonian times a matrix
    Y = zeros(size(X))
    mul!(Y, H, X)
    return Y
    # # new version that acts on the full matrices
    # Y_lap  = lap(H,X)
    # Y_vtot = Vtot(H,X)
    #
    # return Y_lap + Y_vtot
end


function LinearAlgebra.mul!(Y::Array{Float64,2}, H::Hamiltonian, V::Array{Float64,2})
    # in place matrix matrix multiplication
    @assert(size(Y) == size(V))
    Y[:,:] = lap(H,V)
    Y += op_Vtot(H,V)
    return Y
end

# optimized version for eigs (it uses sub arrays to bypass the inference step)
function LinearAlgebra.mul!(
    Y::SubArray{Float64,1,Array{Float64,1}},
    H::Hamiltonian,
    V::SubArray{Float64,1,Array{Float64,1}}
)
    # in place matrix matrix multiplication
    @assert(size(Y) == size(V))
    for ii = 1:size(V,2)
        Y[:,ii] = H*V[:,ii]
    end

end

function size(H::Hamiltonian)
    return (H.Ns, H.Ns)
end

function eltype(H::Hamiltonian)
    # we work always in real space, the Fourier operations are only meant as
    # a pseudo spectral discretization
    return typeof(1.0)
end

function LinearAlgebra.issymmetric(H::Hamiltonian)
    return true
end


include("update_psi.jl")
include("update_rho.jl")

function create_Hamiltonian(H::Hamiltonian)
    # create the matrix version of the Hmailtonian
    A = real(ifft(diagm(0 => H.kmul[:])*fft(Matrix{Float64}(I, length(H.kmul), length(H.kmul)),1),1))
    A += diagm(0 => H.Vtot[:])
    # we symmetrize A
    return 0.5*(A + A')
end

include("op_lap.jl")



function prec(H::Hamiltonian, x::Array{Float64,1})
    # preconditioner for lobpcg
    ## matlab code to translate
    #     for ik = 1:nkpts
    #     X  = gkincell{ik}
    #     Y  = 27.0 + X.*(18.0 + X.*(12.0 + 8.0*X))
    #     p{ik}  = Y./(Y + 16.0*X.^4)
    # end

    # inv_kmul = zeros(size(H.kmul))
    # inv_kmul[1] = 0
    # inv_kmul[2:end] = 1./H.kmul[2:end]
    inv_kmul = zeros(size(H.kmul))
    inv_kmul[1] = 0
    inv_kmul[2:end] = H.kmul[2:end]
    Y = 27.0 .+ inv_kmul.*(18.0 .+ inv_kmul.*(12.0 .+ 8.0 .*inv_kmul))
    yfft = Y./(Y + 16.0 .*inv_kmul.^4)
    ytemp = yfft.*fft(x)
    return real(ifft(ytemp))
end

include("op_Vtot.jl")

include("update_vtot.jl")

import Statistics.mean
using FFTW

function init_pot!(H::Hamiltonian, nocc::Int64)
    # nocc number of occupied states
    #function to initialize the potential in the Hamiltonian class
    rho = -H.rhoa
    rho = rho / ( sum(rho)*H.dx) * (nocc*H.nspin)
    H.rho = rho
    #println("nocc*H.nspin = ", nocc*H.nspin)
    #println("Sum H.rho + H.rhoa = ", sum(H.rho + H.rhoa))
    H.Vhar = hartree_pot_bc(H.rho + H.rhoa, H)
    #println("sum VHar = ", sum(H.Vhar))
    #exit()
    H.Vtot = H.Vhar   # No exchange-correlation
    H.Vtot = H.Vtot .- mean(H.Vtot)# IMPORTANT (zero mean?)
    println("At the end of init_pot: sum(H.Vtot) = ", sum(H.Vtot))
    println("Some H.Vtot:")
    for i in 1:10
        println(H.Vtot[i])
    end
    #exit()
end

function update_pot!(H::Hamiltonian)
     # TODO: I dont' know in how many different ways this si updateded
     # computing the hartree potenatial
    H.Vhar = hartree_pot_bc(H.rho + H.rhoa,H)
    # here Vtotnew only considers the
    Vtotnew  = H.Vhar  # no exchange-correlation so far
    Verr = norm(Vtotnew-H.Vtot)./norm(H.Vtot) # computing the relative error


    # NOTE: H.Fband is only the band energy here.  The real total energy
    # is calculated using the formula below:
    H.Ftot = H.Fband + 1/2 * sum((H.rhoa-H.rho).*H.Vhar)*H.dx
    return (Vtotnew,Verr) # returns the differnece betwen two consecutive iterations
end

include("scf_potmix.jl")


# Function to compute the force and update ham.atoms.force
function get_force!(H::Hamiltonian)
    atoms  = H.atoms
    Natoms = atoms.Natoms
    rhotot = H.rho + H.rhoa
    for i=1:Natoms
        # IMPORTANT: derivative is taken w.r.t atom positions, which introduces the minus sign
        dV = -hartree_pot_bc_opt_vec(H.drhoa[:,i], H.Ls, H.YukawaK, H.epsil0)
        # Force is negative gradient
        atoms.force[i] = - sum(dV.*rhotot)*H.dx
    end
end
