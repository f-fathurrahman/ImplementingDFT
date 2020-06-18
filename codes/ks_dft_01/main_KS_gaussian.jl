push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

function pot_gaussian( grid::FD3dGrid; A=1.0, α=1.0, center=[0.0, 0.0, 0.0], normalized=false )
    Npoints = grid.Npoints
    Vpot = zeros(Npoints)
    if normalized
        NN = 1.0/sqrt(α/pi)^3
    else
        NN = 1.0
    end
    for i in 1:Npoints
        x = grid.r[1,i] - center[1]
        y = grid.r[2,i] - center[2]
        z = grid.r[3,i] - center[3]
        r2 = x^2 + y^2 + z^2
        Vpot[i] = -A*exp( -α*r2 )*NN
    end
    return Vpot
end

function main()
    Random.seed!(1234)

    AA = -3.0*ones(3)
    BB =  3.0*ones(3)
    NN = [25, 25, 25]

    grid = FD3dGrid( NN, AA, BB )

    println("hx = ", grid.hx)
    println("hy = ", grid.hy)
    println("hz = ", grid.hz)
    println("dVol = ", grid.dVol)
    println(grid.hx*grid.hy*grid.hz)

    my_pot_local( grid ) = pot_gaussian( grid, α=1.0, A=1.0, normalized=true )

    Nstates = 1
    Nelectrons = 2*Nstates
    Ham = Hamiltonian( grid, my_pot_local, Nelectrons=Nelectrons )

    Nbasis = prod(NN)

    dVol = grid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi)
    psi = psi/sqrt(dVol)

    for i in 1:Nstates
        @printf("%18.10f\n", dot(psi[:,i], psi[:,i])*dVol )
    end

    Rhoe = zeros(Float64, Nbasis)
    Rhoe_new = zeros(Float64, Nbasis)

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
