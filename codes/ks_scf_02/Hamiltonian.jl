mutable struct Hamiltonian
    grid::Union{FD3dGrid,LF3dGrid}
    Laplacian::SparseMatrixCSC{Float64,Int64}
    V_Ps_loc::Vector{Float64}
    V_Hartree::Vector{Float64}
    V_XC::Vector{Float64}
    pspots::Array{PsPot_GTH}
    pspotNL::PsPotNL
    electrons::Electrons
    atoms::Atoms
    rhoe::Vector{Float64}
    precKin
    precLaplacian
    energies::Energies
end

# FIXME: Not yet used
struct GridInfo
    #
    grid_type::Symbol
    #
    x_min::Float64
    x_max::Float64
    #
    y_min::Float64
    y_max::Float64
    #
    z_min::Float64
    z_max::Float64
    #
    Nx::Int64
    Ny::Int64
    Nz::Int64
end


"""
Build a Hamiltonian with given FD grid and local potential.
"""
function Hamiltonian(
    atoms::Atoms, pspfiles::Array{String,1}, grid;
    Nstates_extra=0,
    func_1d=build_D2_matrix_9pt
)

    # Need better mechanism for this
    @printf("Building Laplacian ...")
    if typeof(grid) == FD3dGrid
        Laplacian = build_nabla2_matrix( grid, func_1d=func_1d )
    else
        Laplacian = build_nabla2_matrix( grid )
    end
    @printf("... done\n")

    @printf("Building preconditioners ...")
    precKin = aspreconditioner( ruge_stuben(-0.5*Laplacian) )
    precLaplacian = aspreconditioner( ruge_stuben(Laplacian) )
    @printf("... done\n")


    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    Npoints = grid.Npoints

    pspots = Array{PsPot_GTH}(undef,Nspecies)
    for isp = 1:Nspecies
        pspots[isp] = PsPot_GTH( pspfiles[isp] )
    end

    V_Ps_loc = zeros(Float64,Npoints)
    atm2species = atoms.atm2species
    for ia in 1:Natoms
        isp = atm2species[ia]
        for ip in 1:Npoints
            dx = grid.r[1,ip] - atoms.positions[1,ia]
            dy = grid.r[2,ip] - atoms.positions[2,ia]
            dz = grid.r[3,ip] - atoms.positions[3,ia]
            dr = sqrt(dx^2 + dy^2 + dz^2)
            V_Ps_loc[ip] = V_Ps_loc[ip] + eval_Vloc_R( pspots[isp], dr )
        end
    end
    println("sum V_Ps_loc = ", sum(V_Ps_loc))

    pspotNL = PsPotNL( atoms, pspots, grid )
    println("NbetaNL = ", pspotNL.NbetaNL)

    V_Hartree = zeros(Float64, Npoints)
    V_XC = zeros(Float64, Npoints)
    rhoe = zeros(Float64, Npoints)

    electrons = Electrons( atoms, pspots, Nstates_empty=Nstates_extra )

    energies = Energies()

    return Hamiltonian( grid, Laplacian, V_Ps_loc, V_Hartree, V_XC, pspots, pspotNL, electrons, atoms,
                        rhoe, precKin, precLaplacian, energies )
end


function op_V_Ps_nloc( Ham::Hamiltonian, psi::Array{Float64,2} )
    Nstates = size(psi,2)

    atoms = Ham.atoms
    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    pspots = Ham.pspots
    prj2beta = Ham.pspotNL.prj2beta
    betaNL = Ham.pspotNL.betaNL

    betaNL_psi = calc_betaNL_psi( Ham.pspotNL.betaNL, psi )
    
    Npoints = Ham.grid.Npoints
    Vpsi = zeros(Float64,Npoints,Nstates)

    for ist = 1:Nstates
        for ia = 1:Natoms
            isp = atm2species[ia]
            psp = pspots[isp]
            for l = 0:psp.lmax
            for m = -l:l
                for iprj = 1:psp.Nproj_l[l+1]
                for jprj = 1:psp.Nproj_l[l+1]
                    ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                    jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                    hij = psp.h[l+1,iprj,jprj]
                    for ip = 1:Npoints
                        Vpsi[ip,ist] = Vpsi[ip,ist] + hij*betaNL[ip,ibeta]*betaNL_psi[ist,jbeta]
                    end
                end # iprj
                end # jprj
            end # m
            end # l
        end
    end
    return Vpsi
end

import Base: *
function *( Ham::Hamiltonian, psi::Matrix{Float64} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    Hpsi = zeros(Float64,Nbasis,Nstates)
    #
    Hpsi = -0.5*Ham.Laplacian * psi
    #
    if Ham.pspotNL.NbetaNL > 0
        Vnlpsi = op_V_Ps_nloc(Ham, psi)
        for ist in 1:Nstates, ip in 1:Nbasis
            Hpsi[ip,ist] = Hpsi[ip,ist] + ( Ham.V_Ps_loc[ip] +
                Ham.V_Hartree[ip] + Ham.V_XC[ip] ) * psi[ip,ist] + Vnlpsi[ip,ist]
        end
    else
        for ist in 1:Nstates, ip in 1:Nbasis
            Hpsi[ip,ist] = Hpsi[ip,ist] + ( Ham.V_Ps_loc[ip] +
                Ham.V_Hartree[ip] + Ham.V_XC[ip] ) * psi[ip,ist]
        end
    end
    return Hpsi
end

function update!( Ham::Hamiltonian, Rhoe::Vector{Float64} )
    Ham.rhoe = Rhoe
    Ham.V_Hartree = Poisson_solve_PCG( Ham.Laplacian, Ham.precLaplacian, -4*pi*Rhoe, 1000, verbose=false, TOL=1e-10 )
    Ham.V_XC = excVWN( Rhoe ) + Rhoe .* excpVWN( Rhoe )
    return
end
