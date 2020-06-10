push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

function pot_harmonic( grid; ω=1.0, center=[0.0, 0.0] )
    Npoints = grid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = grid.r[1,i] - center[1]
        y = grid.r[2,i] - center[2]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 )
    end
    return Vpot
end

function main()

    Random.seed!(1234)

    AA = [-25.0, -25.0]
    BB = [ 25.0,  25.0]
    NN = [81, 81]

    grid = FD2dGrid( NN, AA, BB )

    V_ext = pot_harmonic( grid, ω=0.22 )
    
    Nstates = 1
    Nelectrons = 2*Nstates
    Ham = Hamiltonian( grid, V_ext, Nelectrons=Nelectrons )

    @printf("sizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)

    Nbasis = prod(NN)

    dVol = grid.dVol

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

    #e = 2.0*sum(eigenval(1:N_occ)) + ec + ex
    #e = e - sum(rho(:, :)*(vx(:, :) + vc(:, :) + vh(:, :)/2.0))*delta**2

    println()
    println("calc_energies v2")
    Etot = 2*sum(evals) + sum( excVWN(Ham.rhoe) .* Ham.rhoe )*dVol
    Etot = Etot - sum( Ham.rhoe .* ( 0.5*Ham.V_Hartree .+  Ham.V_XC ) ) * dVol
    println("Etot = ", Etot)
end

@time main()
