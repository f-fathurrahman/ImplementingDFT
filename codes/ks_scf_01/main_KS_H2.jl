push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

include("../common/potential_H_atom.jl")

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
    NN = [50,50,50]

    grid = FD3dGrid( NN, AA, BB )
    println(grid)

    atoms = Atoms(xyz_string="""
    2

    H    0.75   0.0   0.0
    H   -0.75   0.0   0.0
    """)

    V_Ps_loc = pot_Hps_HGH( atoms, grid )

    Nstates = 1
    Nelectrons = 2
    Ham = Hamiltonian( atoms, grid, V_Ps_loc, Nelectrons=Nelectrons )

    Npoints = grid.Npoints
    dVol = grid.dVol

    psi = rand(Float64,Npoints,Nstates)
    ortho_sqrt!(psi)
    psi = psi/sqrt(dVol)

    Rhoe_new = zeros(Float64,Npoints)
    Rhoe = zeros(Float64,Npoints)

    calc_rhoe!( Ham, psi, Rhoe )
    @printf("Integrated Rhoe = %18.10f\n", sum(Rhoe)*dVol)

    update!( Ham, Rhoe )

    evals = zeros(Float64,Nstates)
    Etot_old = 0.0
    dEtot = 0.0
    betamix = 0.5
    dRhoe = 0.0
    NiterMax = 100

    Ham.energies.NN = calc_E_NN( [1.0, 1.0], atoms.positions )

    for iterSCF in 1:NiterMax

        evals = diag_LOBPCG!( Ham, psi, Ham.precKin, verbose_last=false )
        #evals = diag_Emin_PCG!( Ham, psi, Ham.precKin, verbose_last=false )

        #psi = psi*sqrt(dVol) # for diag_davidson
        #evals = diag_davidson!( Ham, psi, Ham.precKin, verbose_last=false )

        psi = psi/sqrt(dVol) # renormalize

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
