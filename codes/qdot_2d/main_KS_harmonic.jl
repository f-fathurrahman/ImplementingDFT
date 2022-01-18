push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

function pot_harmonic( grid; ω=1.0, center=[0.0, 0.0], A=1.0 )
    Npoints = grid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = grid.r[1,i] - center[1]
        y = grid.r[2,i] - center[2]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 ) - A
    end
    return Vpot
end

function main()

    Random.seed!(1234)

    AA = [-15.0, -15.0]
    BB = [ 15.0,  15.0]
    NN = [80, 80]

    grid = FD2dGrid( NN, AA, BB )
    #grid = LF2dGrid( NN, AA, BB, types=(:sinc,:sinc) )

    Npoints = grid.Npoints
    dVol = grid.dVol

    V_ext = pot_harmonic( grid, ω=0.22, A=1.0 )
    
    Nstates = 10
    Nelectrons = 2*Nstates
    Ham = Hamiltonian( grid, V_ext, Nelectrons=Nelectrons )

    @printf("sizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)

    psi = rand(Float64,Npoints,Nstates)
    ortho_sqrt!(psi, dVol)

    Rhoe = zeros(Float64,Npoints)
    Rhoe_new = zeros(Float64,Npoints)

    calc_rhoe!( Ham, psi, Rhoe )
    @printf("Integrated Rhoe = %18.10f\n", sum(Rhoe)*dVol)

    update!( Ham, Rhoe )

    evals = zeros(Float64,Nstates)
    Etot_old = 0.0
    dEtot = 0.0
    betamix = 0.2
    dRhoe = 0.0
    NiterMax = 200
    etot_conv_thr = 1e-6
    Nconverges = 1

    for iterSCF in 1:NiterMax

        @views evals[:] = diag_LOBPCG!( Ham, psi, Ham.precKin, verbose_last=false )
        #evals = diag_Emin_PCG!( Ham, psi, Ham.precKin, verbose_last=false )
        @views psi[:] = psi[:]/sqrt(dVol) # renormalize

        #evals = diag_davidson!( Ham, psi, Ham.precKin, verbose_last=false )
        #psi = psi*sqrt(dVol) # renormalize for diag_davidson

        calc_rhoe!( Ham, psi, Rhoe_new )
        Rhoe = betamix*Rhoe_new + (1-betamix)*Rhoe

        update!( Ham, Rhoe )

        calc_energies!( Ham, psi )
        Etot = sum( Ham.energies )

        dRhoe = sum(abs.(Rhoe - Rhoe_new))/Npoints
        dEtot = abs(Etot - Etot_old)

        @printf("%5d %18.10f %18.10e %18.10e\n", iterSCF, Etot, dEtot, dRhoe)

        if dEtot < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        if Nconverges >= 1
            @printf("\nSCF is converged in iter: %d\n", iterSCF)
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
    Etot = 2*sum(evals) + calc_E_xc_2d( Ham.rhoe, dVol=dVol )
    Etot = Etot - sum( Ham.rhoe .* ( 0.5*Ham.V_Hartree .+  Ham.V_XC ) ) * dVol
    println("Etot = ", Etot)
end

main()
@time main()