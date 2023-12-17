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
    
    println("Initial Haux: ")
    display(Haux); println()

    g = zeros(Npoints,Nstates)
    Kg = zeros(Npoints,Nstates)
    
    Hsub = zeros(Nstates,Nstates)
    g_Haux = zeros(Nstates,Nstates)
    Kg_Haux = zeros(Nstates,Nstates)
    
    d = zeros(Npoints,Nstates)
    d_Haux = zeros(Nstates,Nstates)
    
    
    # Evaluate total energy by calling Lfunc
    E1 = calc_Lfunc_Haux!(Ham, psi, Haux)

    # Evaluate gradients
    calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)
    
    # Precondition
    prec_invK!(Ham, g, Kg)
    
    println("g_Haux = ")
    display(g_Haux)
    #display(Kg_Haux)

    println("dot g,g = ", dot(g,g)*hx)

    dg = dot(g_Haux,g_Haux)*hx
    println("dot g_Haux,g_Haux = ", dg)

    # We use different "learning rate" for psi and Haux
    α = 5e-3 #3e-5
    α_Haux = 0.1

    Urot2 = zeros(Float64, Nstates,Nstates)
    E_new = 0.0
    Udagger = zeros(Float64, Nstates,Nstates)

    for iterEmin in 1:2000
    
        println("----------------")
        println("Enter step: ", iterEmin)
        println("----------------")

        # Set the directions for psi and Haux
        d[:,:] = -Kg[:,:]
        d_Haux[:,:] = -Kg_Haux[:,:]

        # Update wavefunction and auxiliary Hamiltonian
        psi[:,:] = psi[:,:] + α*d[:,:]
        Haux[:,:] = Haux[:,:] + α_Haux*d_Haux[:,:]
    
        # Orthonormalize
        #ortho_sqrt!(psi)
        
        # Orthonormalize wavefunction (involves rotation)
        Udagger[:,:] = inv(sqrt(psi'*psi)) ./ sqrt(hx) # rotation
        psi[:,:] = psi*Udagger

        #println("Check ortho 1:")
        #display(psi' * psi * hx)

        # Also rotate Haux according to Udagger
        Haux[:,:] = Udagger' * Haux * Udagger

        #println("Haux before (after rotated by Udagger):")
        #display(Haux); println()
        
        # This is needed to make Haux diagonal
        # wavefunction psi also must be rotated or transformed
        # Urot2 is the matrix that diagonalizes Haux (probably not needed?)
        Urot2[:,:] = transform_psi_Haux!(psi, Haux)
        # XXX: probably the name of this function should be more specific:
        # We need this because Haux need to be diagonal at all times
        # Calculation of gradient will be wrong if we don't do this
        # The formula used for gradient is assuming that Haux is diagonal
    
        # Evaluate new energy with new psi and Haux
        E_new = calc_Lfunc_Haux!(Ham, psi, Haux)
        # also calculate the gradients
        calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)
        
        prec_invK!(Ham, g, Kg)

        # why????
        #g_Haux[:,:] = Urot2 * g_Haux * Urot2'
        #Kg_Haux[:,:] = Urot2 * Kg_Haux * Urot2'

        #g[:,:] = g*Urot2
        #Kg[:,:] = Kg*Urot2

        println("E1 = ", E1)
        println("E_new = ", E_new)

        # Some heuristics
        if E1 < E_new
            println("WARNING: Energy does not decrease")
            α = α*0.9
            α_Haux = α_Haux*0.9
            println("α and α_Haux are reduced to: ", α, " ", α_Haux)
        end

        dE = abs(E1 - E_new)
        dg_Haux = dot(g_Haux,g_Haux)*hx/length(g_Haux)
        dg = dot(g,g)*hx/length(g)
        println("dot g_Haux,g_Haux = ", dg)

        @printf("Emin: %8d %18.10f %10.5e [%10.5e,%10.5e]\n", iterEmin, E_new, dE, dg, dg_Haux)

        # Check convergence
        # Probably we don't need to set too small criteria for dg
        if (abs(E1 - E_new) < 1e-7) && (dg_Haux < 1e-3)
            println("Converged")
            break
        end

        E1 = E_new # set old value
        Urot[:,:] = Urot2[:,:]

    end

    println("ebands = ")
    display(Ham.electrons.ebands); println()
    
    println("Focc = ")
    println(Ham.electrons.Focc); println()

end
main()