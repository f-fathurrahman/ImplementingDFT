using Printf
using LinearAlgebra
using SparseArrays
using AlgebraicMultigrid
using Random

include("INC_hartree_scf.jl")

function pot_H_atom( atoms::Atoms, grid )
    Npoints = grid.Npoints
    Vpot = zeros(Npoints)
    
    Natoms = atoms.Natoms
    atpos = atoms.positions

    for ia in 1:Natoms, ip in 1:Npoints
        dx = grid.r[1,ip] - atpos[1,ia]
        dy = grid.r[2,ip] - atpos[2,ia]
        dz = grid.r[3,ip] - atpos[3,ia]
        Vpot[ip] = -1.0/sqrt(dx^2 + dy^2 + dz^2)
    end
    return Vpot
end

function pot_Hps_HGH( atoms::Atoms, grid )
    Npoints = grid.Npoints
    Vpot = zeros( Float64, Npoints )

    # Parameters
    Zval = 1
    rloc = 0.2
    C1 = -4.0663326
    C2 = 0.6678322

    Natoms = atoms.Natoms
    atpos = atoms.positions

    for ia in 1:Natoms, ip in 1:Npoints
        dx2 = ( grid.r[1,ip] - atpos[1,ia] )^2
        dy2 = ( grid.r[2,ip] - atpos[2,ia] )^2
        dz2 = ( grid.r[3,ip] - atpos[3,ia] )^2
        r = sqrt(dx2 + dy2 + dz2)
        if r < eps()
            Vpot[ip] = -2*Zval/(sqrt(2*pi)*rloc) + C1
        else
            rrloc = r/rloc
            Vpot[ip] = -Zval/r * erf( r/(sqrt(2.0)*rloc) ) +
                     (C1 + C2*rrloc^2)*exp(-0.5*(rrloc)^2)
        end
    end
    return Vpot
end

function main()

    Random.seed!(1234)

    AA = [-8.0, -8.0, -8.0]
    BB = [ 8.0,  8.0,  8.0]
    NN = [50, 50, 50]

    grid = FD3dGrid( NN, AA, BB )

    atoms = Atoms( xyz_string=
        """
        1

        H  0.0  0.0  0.0
        """ )

    V_Ps_loc = pot_Hps_HGH(atoms, grid)

    Nstates = 1
    Nelectrons = 1

    Ham = Hamiltonian(atoms, grid, V_Ps_loc, Nelectrons=1)

    Npoints = grid.Npoints
    dVol = grid.dVol

    psi = rand(Float64,Npoints,Nstates)
    ortho_sqrt!(psi)
    psi = psi/sqrt(dVol)

    println("Test normalization")
    for i in 1:Nstates
        @printf("%18.10f\n", dot(psi[:,i], psi[:,i])*dVol )
    end

    Rhoe = calc_rhoe( Ham, psi )
    @printf("Integrated Rhoe = %18.10f\n", sum(Rhoe)*dVol)

    update!( Ham, Rhoe )

    evals = zeros(Float64,Nstates)
    Etot_old = 0.0
    dEtot = 0.0
    betamix = 0.5
    dRhoe = 0.0
    NiterMax = 100

    for iterSCF in 1:NiterMax

        evals = diag_LOBPCG!( Ham, psi, Ham.precKin, verbose_last=false )
        #evals = diag_Emin_PCG!( Ham, psi, Ham.precKin, verbose_last=true )
        psi = psi/sqrt(dVol)

        Rhoe_new = calc_rhoe( Ham, psi )
        #@printf("Integrated Rhoe_new = %18.10f\n", sum(Rhoe_new)*dVol)

        Rhoe = betamix*Rhoe_new + (1-betamix)*Rhoe
        #@printf("Integrated Rhoe     = %18.10f\n", sum(Rhoe)*dVol)

        update!( Ham, Rhoe )

        Etot = sum( calc_energies( Ham, psi ) )

        dRhoe = sum(abs.(Rhoe - Rhoe_new))/Npoints # MAE
        dEtot = abs(Etot - Etot_old)

        @printf("%5d %18.10f %18.10e %18.10e\n", iterSCF, Etot, dEtot, dRhoe)

        if dEtot < 1e-6
            @printf("Convergence is achieved in %d iterations\n", iterSCF)
            for i in 1:Nstates
                @printf("%3d %18.10f\n", i, evals[i])
            end
            break
        end

        Etot_old = Etot
    end
end

main()
