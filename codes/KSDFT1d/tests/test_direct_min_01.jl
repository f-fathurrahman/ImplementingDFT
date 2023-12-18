push!(LOAD_PATH, "../")

using Printf
using LinearAlgebra
using KSDFT1d
import Random

include("system_defs_01.jl")
include("../utilities.jl")
include("direct_min_no_smearing.jl")


function main()

    Random.seed!(1234)

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

    for iterEmin in 1:50
        
        Etot_old = Etot
        Etot = calc_KohnSham_Etotal!(Ham, psi)

        g[:,:] .= calc_grad(Ham, psi)

        dg = dot(g,g)*hx
        dEtot = abs(Etot - Etot_old)
        @printf("iterEmin = %3d Etot = %18.10f dEtot = %10.5e dg = %10.5e\n", 
            iterEmin, Etot, dEtot, dg)
        if dEtot < 1e-6
            println("Converged")
            break
        end

        Kg[:,:] .= prec_invK(Ham, g)
        #Kg[:,:] .= prec_invHam(Ham, g)

        if iterEmin >= 2
            β = real( dot(g - g_old, Kg) )/real( dot(g_old,Kg_old) )
            if β < 0.0
                β = 0.0
            end
            #println("β = ", β)
        end

        @views d[:] .= -Kg[:] + β*d_old[:]
        psic[:] = ortho_sqrt( psi + α_t*d )
        psic[:] = psic[:]/sqrt(hx)

        #_ = calc_KohnSham_Etotal!(Ham, psic)
        gt[:,:] = calc_grad(Ham, psic; update=true)

        # conj is not really needed because wavefunctions are real-valued here
        denum = real(sum(conj(g-gt).*d))
        if denum != 0.0
            α = abs( α_t*real(sum(conj(g).*d))/denum )
        else
            α = 0.0
        end

        # Update wavefunction
        psi[:,:] .= ortho_sqrt(psi + α*d)
        psi[:] = psi[:]/sqrt(hx)

        # Previous
        @views g_old[:] .= g[:]
        @views Kg_old[:] .= Kg[:]
        @views d_old[:] .= d[:]
    end

    Hr = psi' * (Ham*psi) * hx
    band_energies = eigvals(Hr)
    display(band_energies); println()

    println(Ham.energies)

end

main()
#@time main()