push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

function main()
    
    atoms = Atoms( xyz_file="NH3.xyz" )
    println(atoms)

    pspfiles = ["N-q5.gth", "H-q1.gth"]

    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [41, 41, 41]

    grid = FD3dGrid( NN, AA, BB )
    #grid = LF3dGrid( NN, AA, BB )
    
    println("hx = ", grid.hx)
    println("hy = ", grid.hy)
    println("hz = ", grid.hz)
    println("dVol = ", grid.dVol)
    println(grid.hx*grid.hy*grid.hz)

    Ham = Hamiltonian( atoms, pspfiles, grid )

    Nbasis = prod(NN)
    Nstates = Ham.electrons.Nstates
    dVol = grid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi)
    psi = psi/sqrt(dVol)

    for i in 1:Nstates
        @printf("%18.10f\n", dot(psi[:,i], psi[:,i])*dVol )
    end

    Rhoe_new = zeros(Float64,Nbasis)
    Rhoe = zeros(Float64,Nbasis)
    
    calc_rhoe!( Ham, psi, Rhoe )
    @printf("Integrated Rhoe = %18.10f\n", sum(Rhoe)*dVol)

    println(Ham.electrons)

    update!( Ham, Rhoe )
    
    evals = zeros(Float64,Nstates)
    Etot_old = 0.0
    dEtot = 0.0
    betamix = 0.5
    dRhoe = 0.0
    NiterMax = 100

    ethr_evals_last=1e-5
    ethr = 0.1

    Ham.energies.NN = calc_E_NN( Ham.atoms, Ham.pspots )

    for iterSCF in 1:NiterMax

        # determine convergence criteria for diagonalization
        if iterSCF == 1
            ethr = 0.1
        elseif iterSCF == 2
            ethr = 0.01
        else
            ethr = ethr/5.0
            ethr = max( ethr, ethr_evals_last )
        end

        evals = diag_LOBPCG!( Ham, psi, Ham.precKin, verbose_last=true )
        #evals = diag_LOBPCG!( Ham, psi, Ham.precKin, verbose=true, tol=ethr )
        #evals = diag_Emin_PCG!( Ham, psi, Ham.precKin, verbose_last=false )

        #psi = psi*sqrt(dVol) # for diag_davidson
        #evals = diag_davidson!( Ham, psi, Ham.precKin, verbose_last=false )

        psi = psi/sqrt(dVol) # renormalize

        #Rhoe_new = calc_rhoe( Ham, psi )
        calc_rhoe!( Ham, psi, Rhoe_new )
        @printf("Integrated Rhoe              = %18.10f\n", sum(Rhoe)*dVol)

        Rhoe = betamix*Rhoe_new + (1-betamix)*Rhoe

        @printf("Integrated Rhoe after mixing = %18.10f\n", sum(Rhoe)*dVol)

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

main()