push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

function pot_Hps_HGH( grid, center )
    Npoints = grid.Npoints
    Vpot = zeros( Float64, Npoints )

    # Parameters
    Zval = 1
    rloc = 0.2
    C1 = -4.0663326
    C2 = 0.6678322

    # TODO Add journal reference
    for ip = 1:Npoints
        r = norm( grid.r[:,ip] - center[:] )
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

function calc_E_NN( Zvals::Array{Float64,1}, r::Array{Float64,2} )
    Natoms = length(Zvals)
    @assert Natoms == size(r,2)

    E_NN = 0.0
    for ia in 1:Natoms
        for ja in ia+1:Natoms
            dx = r[1,ja] - r[1,ia]
            dy = r[2,ja] - r[2,ia]
            dz = r[3,ja] - r[3,ia]
            r_ij = sqrt(dx^2 + dy^2 + dz^2)
            E_NN = E_NN + Zvals[ia]*Zvals[ja]/r_ij
        end
    end

    return E_NN
end

function main()
    Random.seed!(1234)

    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [65, 65, 65]

    grid = FD3dGrid( NN, AA, BB )

    println("hx = ", grid.hx)
    println("hy = ", grid.hy)
    println("hz = ", grid.hz)
    println("dVol = ", grid.dVol)
    println(grid.hx*grid.hy*grid.hz)


    r_1 = [ 0.75, 1e-9, 1e-9]
    Vion_1 = pot_Hps_HGH(grid, r_1)
    
    r_2 = [-0.75, 1e-9, 1e-9]
    Vion_2 = pot_Hps_HGH(grid, r_2)
    
    my_pot_local = Vion_1 .+ Vion_2

    Nstates = 1
    Nelectrons = 2
    Ham = Hamiltonian( grid, my_pot_local, Nelectrons=Nelectrons, func_1d=build_D2_matrix_9pt )

    Nbasis = prod(NN)

    dVol = grid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi)
    psi = psi/sqrt(dVol)

    for i in 1:Nstates
        @printf("%18.10f\n", dot(psi[:,i], psi[:,i])*dVol )
    end

    Rhoe_new = zeros(Float64,Nbasis)
    Rhoe = zeros(Float64,Nbasis)

    #Rhoe = calc_rhoe( Ham, psi )
    calc_rhoe!( Ham, psi, Rhoe )
    @printf("Integrated Rhoe = %18.10f\n", sum(Rhoe)*dVol)

    update!( Ham, Rhoe )

    evals = zeros(Float64,Nstates)
    Etot_old = 0.0
    dEtot = 0.0
    betamix = 0.5
    dRhoe = 0.0
    NiterMax = 100

    Natoms = 2
    r = zeros(3,Natoms)
    r[:,1] = r_1[:]
    r[:,2] = r_2[:]
    Ham.energies.NN = calc_E_NN( [1.0, 1.0], r )

    for iterSCF in 1:NiterMax

        evals = diag_LOBPCG!( Ham, psi, Ham.precKin, verbose_last=false )
        #evals = diag_Emin_PCG!( Ham, psi, Ham.precKin, verbose_last=false )

        #psi = psi*sqrt(dVol) # for diag_davidson
        #evals = diag_davidson!( Ham, psi, Ham.precKin, verbose_last=false )

        psi = psi/sqrt(dVol) # renormalize

        #Rhoe_new = calc_rhoe( Ham, psi )
        calc_rhoe!( Ham, psi, Rhoe_new )

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
