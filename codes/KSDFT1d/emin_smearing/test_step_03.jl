push!(LOAD_PATH, "../")

import Random
using Printf
using LinearAlgebra
using Serialization

using KSDFT1d

include("system_defs_01.jl")

include("Lfunc.jl")
include("../utilities.jl")
include("gradients_psi_Haux.jl")


# TODO: move out Hamiltonian update steps from calc_Lfunc! and calc_grad_Lfunc!

function main()

    # Initialize a Hamiltonian object
    Ham = init_Hamiltonian()
    
    hx = Ham.grid.hx
    Npoints = Ham.grid.Npoints
    Nelectrons = Ham.electrons.Nelectrons
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc

    # Random wavefunc    
    Random.seed!(1234)
    psi = generate_random_wavefunc(Ham)

    update_from_wavefunc!(Ham, psi) # update the potential
    
    # Prepare Haux from psi (form subspace Hamiltonian)
    #Haux = psi' * (Ham*psi) * hx # Hsub, subspace Hamiltonian
    
    # Using random diagonal Haux
    ebands1 = sort(randn(Nstates))
    Haux = diagm(0 => ebands1)
    
    # Only needed if Haux is not diagonal
    Urot = transform_psi_Haux!(psi, Haux)

    g = zeros(Npoints,Nstates)
    Kg = zeros(Npoints,Nstates)
    
    Hsub = zeros(Nstates,Nstates)
    g_Haux = zeros(Nstates,Nstates)
    Kg_Haux = zeros(Nstates,Nstates)
    
    d = zeros(Npoints,Nstates)
    d_Haux = zeros(Nstates,Nstates)

    Urot2 = zeros(Float64, Nstates,Nstates)
    E1 = 1e10 # quite large number
    E_new = 0.0
    Udagger = zeros(Float64, Nstates,Nstates)


    psic  = zeros(Float64, Npoints, Nstates)
    Hauxc = zeros(Float64, Nstates, Nstates)
    gt = zeros(Float64, Npoints, Nstates)
    gt_Haux = zeros(Float64, Nstates, Nstates)
    Kgt_Haux = zeros(Float64, Nstates, Nstates)    

    g_old = zeros(Float64, Npoints, Nstates)
    Kg_old = zeros(Float64, Npoints, Nstates)
    d_old = zeros(Float64, Npoints, Nstates)

    g_Haux_old = zeros(Float64, Nstates, Nstates)
    Kg_Haux_old = zeros(Float64, Nstates, Nstates)
    d_Haux_old = zeros(Float64, Nstates, Nstates)

    α_t = 1e-1
    α_t_Haux = 1e-1
    α = 1e-2
    α_Haux = 1e-1
    β = 0.0
    β_Haux = 0.0

    for iterEmin in 1:100

        println("\niterEmin = $(iterEmin)")

        # Evaluate total energy by calling Lfunc
        E_new = calc_Lfunc_Haux!(Ham, psi, Haux)
        #
        # Evaluate gradients
        calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)
        #
        # Precondition the gradient for psi
        prec_invK!(Ham, g, Kg)
        # Preconditioning for g_Haux is done in calc_grad_Lfunc_Haux!

        dE = abs(E1 - E_new)
        dg = dot(g,g)*hx /  length(g)
        dg_Haux = dot(g_Haux,g_Haux)*hx / length(g_Haux)
        #
        @printf("Emin: %8d %18.10f %10.5e [%10.5e %10.5e]\n", iterEmin, E_new, dE, dg, dg_Haux)

        # Check convergence
        # Probably we don't need to set too small criteria for dg_Haux
        if (abs(E1 - E_new) < 1e-7) && (dg_Haux < 1e-3)
            println("Converged")
            break
        end

        if (iterEmin >= 10) && ( ( abs(α_t) + abs(α_t_Haux) ) < 10*eps() )
            println("STOPPED: too small α")
            break
        end

        # No convergence yet if we reach this point

        # Some heuristics for setting α and α_Haux
        if E1 < E_new
            println("WARNING: Energy does not decrease")
            # activate this heuristic for steepest descent without line minimization
            #α = α*0.9
            #α_Haux = α_Haux*0.9
            #println("α and α_Haux are reduced to: $(α), $(α_Haux)")
            #
            α_t = α_t*0.5
            α_t_Haux = α_t_Haux*0.5
            println("α_t and α_t_Haux are reduced to: $(α_t), $(α_t_Haux)")
        end

#=
        if iterEmin >= 2
            β = real( dot(g - g_old, Kg) )/real( dot(g_old,Kg_old) )
            if β < 0.0
                β = 0.0
            end
            β_Haux = real( dot(g_Haux - g_Haux_old, Kg_Haux) )/real( dot(g_Haux_old,Kg_Haux_old) )
            if β_Haux < 0.0
                β_Haux = 0.0
            end
        end
=#

        println("β = $(β) β_Haux = $(β_Haux)")


        #
        # Set the directions for psi and Haux
        #
        d[:,:] = -Kg[:,:] + β*d_old[:,:]
        d_Haux[:,:] = -Kg_Haux[:,:] + β_Haux*d_Haux_old[:,:]


        #
        # Line mininization to find α
        #
        psic[:,:] = psi[:,:] + α_t*d[:,:] # trial wavefunc
        Hauxc[:,:] = Haux[:,:] + α_t_Haux*d_Haux[:,:]
        #
        Udagger[:,:] = inv(sqrt(psic'*psic)) ./ sqrt(hx) # rotation
        psic[:,:] = psic*Udagger # orthogonalize
        Hauxc[:,:] = Udagger' * Hauxc * Udagger # rotate Haux
        Urot2[:,:] = transform_psi_Haux!(psic, Hauxc) # make Haux diagonal 
        #display(Hauxc)
        #
        _ = calc_Lfunc_Haux!(Ham, psic, Hauxc)
        calc_grad_Lfunc_Haux!(Ham, psic, Hauxc, gt, Hsub, gt_Haux, Kgt_Haux)
        #
        denum = real(sum(conj(g-gt).*d)) # no need for factor of hx for dot product here
        println("denum = ", denum)
        if denum != 0.0
            α = abs( α_t*real(sum(conj(g).*d))/denum )
        else
            α = 0.0
        end
        #
        denum = real(sum(conj(g_Haux-gt_Haux).*d_Haux))
        println("denum = ", denum)
        if denum != 0.0
            α_Haux = abs( α_t_Haux*real(sum(conj(g_Haux).*d_Haux))/denum )
        else
            α_Haux = 0.0
        end
        println("linmin: α = $(α)  α_Haux=$(α_Haux)")


        #
        # Update wavefunction and auxiliary Hamiltonian
        #
        psi[:,:] = psi[:,:] + α*d[:,:]
        Haux[:,:] = Haux[:,:] + α_Haux*d_Haux[:,:]
        #
        # Orthonormalize wavefunction (involves rotation)
        Udagger[:,:] = inv(sqrt(psi'*psi)) ./ sqrt(hx) # rotation
        psi[:,:] = psi*Udagger
        #
        # Also rotate Haux according to Udagger
        Haux[:,:] = Udagger' * Haux * Udagger
        #
        # This is needed to make Haux diagonal
        # wavefunction psi also must be rotated or transformed
        # Urot2 is the matrix that diagonalizes Haux (probably not needed?)
        Urot2[:,:] = transform_psi_Haux!(psi, Haux)
        # XXX: probably the name of this function should be more specific:

        # set old values before going to new iteration
        E1 = E_new
        #Urot[:,:] = Urot2[:,:] # XXXX Need this?
        @views g_old[:] .= g[:]
        @views Kg_old[:] .= Kg[:]
        @views d_old[:] .= d[:]

        @views g_Haux_old[:] .= g_Haux[:]
        @views Kg_Haux_old[:] .= Kg_Haux[:]
        @views d_Haux_old[:] .= d_Haux[:]

    end

    println("ebands = ")
    display(Ham.electrons.ebands); println()
    
    println("Focc = ")
    println(Ham.electrons.Focc); println()

end
main()