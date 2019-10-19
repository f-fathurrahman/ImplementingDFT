push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

function pot_harmonic( fdgrid::FD3dGrid; ω=1.0, center=[0.0, 0.0, 0.0] )
    Npoints = fdgrid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = fdgrid.r[1,i] - center[1]
        y = fdgrid.r[2,i] - center[2]
        z = fdgrid.r[3,i] - center[3]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 + z^2 )
    end
    return Vpot
end

function main()
    Random.seed!(1234)

    AA = [-3.0, -3.0, -3.0]
    BB = [3.0, 3.0, 3.0]
    NN = [25, 25, 25]

    fdgrid = FD3dGrid( NN, AA, BB )

    println("hx = ", fdgrid.hx)
    println("hy = ", fdgrid.hy)
    println("hz = ", fdgrid.hz)
    println("dVol = ", fdgrid.dVol)
    println(fdgrid.hx*fdgrid.hy*fdgrid.hz)

    my_pot_harmonic( fdgrid ) = pot_harmonic( fdgrid, ω=2 )

    Nstates = 4
    Nelectrons = 2*Nstates
    Ham = Hamiltonian( fdgrid, my_pot_harmonic, Nelectrons=Nelectrons, func_1d=build_D2_matrix_9pt )

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
    betamix = 0.8
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
