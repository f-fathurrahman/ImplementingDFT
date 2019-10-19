using Printf
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using IncompleteLU
using AlgebraicMultigrid
using Random
using SpecialFunctions

include("../../3d/FD3dGrid.jl")
include("../../3d/build_nabla2_matrix.jl")
include("../../diag_Emin_PCG.jl")
include("../../diag_davidson.jl")
include("../../diag_LOBPCG.jl")
include("../../ortho_sqrt.jl")
include("../../supporting_functions.jl")
include("../../3d_poisson/Poisson_solve_PCG.jl")

include("Electrons.jl")
include("Energies.jl")
include("Hamiltonian.jl")

include("calc_rhoe.jl")

include("../LDA_VWN.jl")

include("calc_energies.jl")

function pot_Hps_HGH( fdgrid, center )
    Npoints = fdgrid.Npoints
    Vpot = zeros( Float64, Npoints )

    # Parameters
    Zval = 1
    rloc = 0.2
    C1 = -4.0663326
    C2 = 0.6678322

    # TODO Add journal reference
    for ip = 1:Npoints
        r = norm( fdgrid.r[:,ip] - center[:] )
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

    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [50, 50, 50]

    fdgrid = FD3dGrid( NN, AA, BB )

    println("hx = ", fdgrid.hx)
    println("hy = ", fdgrid.hy)
    println("hz = ", fdgrid.hz)
    println("dVol = ", fdgrid.dVol)
    println(fdgrid.hx*fdgrid.hy*fdgrid.hz)

    my_pot_local( fdgrid ) = pot_Hps_HGH(fdgrid, zeros(3))

    Nstates = 1
    Nelectrons = 1
    Ham = Hamiltonian( fdgrid, my_pot_local, Nelectrons=1, func_1d=build_D2_matrix_9pt )

    Nbasis = prod(NN)

    dVol = fdgrid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi)
    psi = psi/sqrt(dVol)

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

        #psi = psi*sqrt(dVol) # for diag_davidson
        #evals = diag_davidson!( Ham, psi, Ham.precKin, verbose_last=false )

        psi = psi/sqrt(dVol) # renormalize

        Rhoe_new = calc_rhoe( Ham, psi )

        Rhoe = betamix*Rhoe_new + (1-betamix)*Rhoe

        update!( Ham, Rhoe )

        calc_energies!( Ham, psi )
        Etot = sum( Ham.energies )

        dRhoe = norm(Rhoe - Rhoe_new)
        dEtot = abs(Etot - Etot_old)

        @printf("%5d %18.10f %18.10e %18.10e\n", iterSCF, Etot, dEtot, dRhoe)

        if dEtot < 1e-6
            @printf("Convergence is achieved in %d iterations\n", iterSCF)
            @printf("\nEigenvalues:\n")
            for i in 1:Nstates
                @printf("%3d %18.10f\n", i, evals[i])
            end
            break
        end

        Etot_old = Etot
    end

    println(Ham.energies)

end

@time main()
@time main()
