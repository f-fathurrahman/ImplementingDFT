const AMG_PREC_TYPE = typeof( aspreconditioner(ruge_stuben(speye(1))) )

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
    precKin::Union{AMG_PREC_TYPE,ILU0Preconditioner}
    psolver::Union{PoissonSolverDAGE,PoissonSolverFFT}
    energies::Energies
    gvec::Union{Nothing,GVectors}
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

function init_V_Ps_loc( atoms::Atoms, grid, pspots::Array{PsPot_GTH,1} )

    @assert atoms.pbc == (false,false,false)

    Npoints = grid.Npoints
    V_Ps_loc = zeros(Float64,Npoints)
    
    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    
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
    return V_Ps_loc
end



function init_V_Ps_loc_G( atoms, grid, gvec, pspots )

    @assert atoms.pbc == (true,true,true)

    Npoints = grid.Npoints
    CellVolume = grid.Lx * grid.Ly * grid.Lz  # FIXME: orthogonal LatVecs
    Nx = grid.Nx
    Ny = grid.Ny
    Nz = grid.Nz
    atm2species = atoms.atm2species
    Nspecies = atoms.Nspecies
    G2 = gvec.G2
    Ng = length(G2)

    V_Ps_loc = zeros(Float64, Npoints)
    Vg = zeros(ComplexF64, Nx,Ny,Nz)

    if typeof(grid) == LF3dGrid
        # periodic LF
        shifts = [0.5*grid.hx, 0.5*grid.hy, 0.5*grid.hz]
        strf = calc_strfact_shifted( atoms, gvec, shifts )
    else
        # periodic FD    
        strf = calc_strfact( atoms, gvec )
    end

    for isp = 1:Nspecies
        psp = pspots[isp]
        for ig = 1:Ng
            Vg[ig] = strf[ig,isp] * eval_Vloc_G( psp, G2[ig] )
        end
        #
        ifft!(Vg)
        @views V_Ps_loc[:] = V_Ps_loc[:] + real( Vg[:] ) * Npoints / CellVolume
    end

    return V_Ps_loc
end


"""
Build a Hamiltonian with given FD grid and local potential.
"""
function Hamiltonian(
    atoms::Atoms, pspfiles::Array{String,1}, grid;
    Nstates_extra=0,
    verbose=true,
    stencil_order=9,
    prec_type=:ILU0
)

    # Need better mechanism for this
    verbose && @printf("Building Laplacian ...")
    if typeof(grid) == FD3dGrid
        Laplacian = build_nabla2_matrix( grid, stencil_order=stencil_order )
    else
        Laplacian = build_nabla2_matrix( grid )
    end
    verbose && @printf("... done\n")

    # Initialize gvec for periodic case
    if atoms.pbc == (true,true,true)
        gvec = GVectors(grid)
    else
        gvec = nothing
    end

    verbose && @printf("Building preconditioners ...")
    if prec_type == :amg
        precKin = aspreconditioner( ruge_stuben(-0.5*Laplacian) )
    else
        precKin = ILU0Preconditioner(-0.5*Laplacian)
    end
    verbose && @printf("... done\n")

    Nspecies = atoms.Nspecies
    Npoints = grid.Npoints

    # Pseudopotentials
    pspots = Array{PsPot_GTH}(undef,Nspecies)
    for isp = 1:Nspecies
        pspots[isp] = PsPot_GTH( pspfiles[isp] )
    end

    # Local pseudopotential
    if atoms.pbc == (true,true,true)
        V_Ps_loc = init_V_Ps_loc_G(atoms, grid, gvec, pspots)
    else
        V_Ps_loc = init_V_Ps_loc(atoms, grid, pspots)
    end

    verbose && println("sum V_Ps_loc = ", sum(V_Ps_loc))
    pspotNL = PsPotNL( atoms, pspots, grid )

    V_Hartree = zeros(Float64, Npoints)
    V_XC = zeros(Float64, Npoints)
    rhoe = zeros(Float64, Npoints)

    electrons = Electrons( atoms, pspots, Nstates_empty=Nstates_extra )

    energies = Energies()

    if atoms.pbc == (false,false,false)
        psolver = PoissonSolverDAGE(grid)
    else
        psolver = PoissonSolverFFT(grid)
    end

    return Hamiltonian( grid, Laplacian, V_Ps_loc, V_Hartree, V_XC, pspots, pspotNL, electrons, atoms,
                        rhoe, precKin, psolver, energies, gvec )
end


function op_V_Ps_nloc( Ham::Hamiltonian, psi::Array{Float64,2} )
    Nstates = size(psi,2)

    atoms = Ham.atoms
    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    pspots = Ham.pspots
    prj2beta = Ham.pspotNL.prj2beta
    betaNL = Ham.pspotNL.betaNL

    dVol = Ham.grid.dVol
    betaNL_psi = psi' * Ham.pspotNL.betaNL *dVol
    
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

function op_H(Ham, psi)
    return Ham*psi
end

function update!( Ham::Hamiltonian, Rhoe::Vector{Float64} )
    Ham.rhoe[:] = Rhoe[:]
    Ham.V_Hartree = Poisson_solve( Ham.psolver, Ham.grid, Rhoe )
    Ham.V_XC = excVWN( Rhoe ) + Rhoe .* excpVWN( Rhoe )
    return
end
