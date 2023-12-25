push!(LOAD_PATH, "../")

using Printf
using LinearAlgebra
using KSDFT1d
import Random

include("../tests/system_defs_01.jl")
include("../utilities.jl")
include("direct_min_no_smearing.jl")


function main()

    #iseed = abs(rand(Int64))
    iseed = 1234
    Random.seed!(iseed)
    println("Using iseed = ", iseed)

    Ham = init_Hamiltonian()

    hx = Ham.grid.hx
    Npoints = Ham.grid.Npoints
    Nelectrons = Ham.electrons.Nelectrons
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc

    Kmat = Ham.Kmat
    Vtot = Ham.potentials.Total
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    rhoe = Ham.rhoe

    Nspin = 1
    Hmat = zeros(Float64, Npoints, Npoints)
    epsxc = zeros(Float64, Npoints)
    rhoe_new = zeros(Float64, Npoints, Nspin)

    Etot = Inf
    Etot_old = Etot
    E_NN = calc_E_NN(Ham.atoms)
    Ham.energies.NN = E_NN

    psi = ortho_sqrt( rand(Float64, Npoints, Nstates) )
    psi[:] = psi[:]/sqrt(hx)
    display(psi' * psi * hx)


    α_t = 3e-5
    g = zeros(Float64, Npoints, Nstates)
    Kg = zeros(Float64, Npoints, Nstates)
    gt = zeros(Float64, Npoints, Nstates)
    d = zeros(Float64, Npoints, Nstates)
    psic = zeros(Float64, Npoints, Nstates)

    g_old = zeros(Float64, Npoints, Nstates)
    Kg_old = zeros(Float64, Npoints, Nstates)
    d_old = zeros(Float64, Npoints, Nstates)

    β = 0.0
    Nconverged = 0

    for iterEmin in 1:100
        
        println("\nBegin iterEmin = ", iterEmin)

        Etot_old = Etot
        # Hamiltonian will be updated by calling calc_KohnSham_Etotal!
        Etot = calc_KohnSham_Etotal!(Ham, psi)
        # Then, we call gradient without updating Hamiltonian
        g[:,:] .= calc_grad(Ham, psi)

        dg = 2*dot(g,g)*hx
        dEtot = abs(Etot - Etot_old)
        @printf("%3d Etot=%18.10f dE=%10.5e dG=%10.5e\n", iterEmin, Etot, dEtot, dg)
        #
        if dEtot < 1e-6
            Nconverged += 1
        else
            Nconverged = 0
        end
        if Nconverged >= 2
            println("Converged")
            break
        end

        prec_invK!(Ham, g, Kg)
        #prec_invHam!(Ham, g, Kg)

        num1 = 2*dot(g,d_old)*hx
        gg = 2*dot(g,g)*hx
        dd = 2*dot(d_old,d_old)*hx
        println("num1 = ", num1)
        println("gg*dd = ", gg*dd)
        linmin_test = num1/sqrt(gg*dd)
        println("linmin_test = ", linmin_test)


        if iterEmin >= 2
            β = dot(g-g_old, Kg)/dot(g_old,Kg_old)
            if β < 0.0
                β = 0.0
            end
        end
        println("β = ", β)

        # new gradient is already evaluated?
        num1 = 2*dot(g, Kg_old)*hx
        gg1 = 2*dot(g, Kg)*hx
        gg2 = 2*dot(g_old, Kg_old)*hx
        println("num1 = ", num1, " gg1 = ", gg1, " gg2 = ", gg2)
        cg_test = num1/sqrt(gg1*gg2)
        println("cg_test = ", cg_test)

        # Set search direction
        d[:] .= -Kg[:] + β*d_old[:]
        constrain_search_dir!(d, psi, hx)
        
        # Line minimization
        psic[:] = ortho_sqrt( psi + α_t*d )
        psic[:] = psic[:]/sqrt(hx)
        #_ = calc_KohnSham_Etotal!(Ham, psic)
        gt[:,:] = calc_grad(Ham, psic; update=true)
        #
        denum = 2*dot(g-gt, d)*hx
        if denum != 0.0
            num = 2*dot(g, d)*hx
            α = abs( α_t*num/denum )
        else
            α = 0.0
        end
        println("α = ", α)

        # Update wavefunction
        psi[:,:] .= ortho_sqrt(psi + α*d)
        psi[:] = psi[:]/sqrt(hx)

        # Previous
        @views g_old[:] .= g[:]
        @views Kg_old[:] .= Kg[:]
        @views d_old[:] .= d[:]
    end

    Hr = psi' * (Ham*psi) * hx
    band_energies, U = eigen(Hr)
    psi[:,:] = psi*U

    println("Band energies:")
    display(band_energies); println()

    println("Total energy components")
    println(Ham.energies)

end

main()
#@time main()