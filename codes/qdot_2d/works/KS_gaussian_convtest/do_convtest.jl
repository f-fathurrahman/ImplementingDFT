# Only tested on Linux
push!(LOAD_PATH, "/home/efefer/WORKS/my_github_repos/ImplementingDFT/codes/qdot_2d/")

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

function pot_gaussian( grid::FD2dGrid; A=1.0, α=1.0 )
    Npoints = grid.Npoints
    f = zeros(Npoints)
    for i in 1:Npoints
        x = grid.r[1,i]
        y = grid.r[2,i]
        r2 = x^2 + y^2
        f[i] = -A*exp(-α*r2)
    end
    return f
end

function do_run(N::Int64; Nstates=1, basis=:FD)

    Random.seed!(1234)

    AA = [-8.0, -8.0]
    BB = [ 8.0,  8.0]
    NN = [N, N]

    if basis == :FD
        grid = FD2dGrid( NN, AA, BB )
    elseif basis == :sinc
        grid = LF2dGrid( NN, AA, BB, types=(:sinc,:sinc) )
    elseif basis == :pib
        grid = LF2dGrid( NN, AA, BB, types=(:C,:C) )
    else
        error("Unknown basis: ", basis)
    end

    V_ext = pot_gaussian( grid, α=1.0, A=10.0 )
    
    Nelectrons = 2*Nstates
    Ham = Hamiltonian( grid, V_ext, Nelectrons=Nelectrons )
    @printf("sizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)

    Npoints = grid.Npoints
    dVol = grid.dVol
    psi = rand(Float64,Npoints,Nstates)
    ortho_sqrt!(psi)
    psi = psi/sqrt(dVol)

    Rhoe = calc_rhoe( Ham, psi )
    @printf("Integrated Rhoe = %18.10f\n", sum(Rhoe)*dVol)

    update!( Ham, Rhoe )

    evals = zeros(Float64,Nstates)
    Etot_old = 0.0
    dEtot = 0.0
    betamix = 0.1
    dRhoe = 0.0
    NiterMax = 200
    etot_conv_thr = 1e-6
    Nconverges = 1

    for iterSCF in 1:NiterMax

        evals = diag_LOBPCG!( Ham, psi, Ham.precKin, verbose_last=false )
        psi = psi/sqrt(dVol) # renormalize

        Rhoe_new = calc_rhoe( Ham, psi )

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
            break
        end

        Etot_old = Etot
    end

    @printf("\nEigenvalues:\n")
    for i in 1:Nstates
        @printf("%3d %18.10f\n", i, evals[i])
    end

    println(Ham.energies)

end

do_run(30, Nstates=2, basis=:FD)
#for N in range(20, stop=80, step=10)
#    do_run(N, Nstates=3, basis=:FD)
#end