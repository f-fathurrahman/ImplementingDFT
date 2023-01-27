push!(LOAD_PATH, "../")

import Random
using Printf
using LinearAlgebra
using Serialization

using KSDFT1d

#include("system_defs_01.jl")
include("system_defs_02.jl")

include("Lfunc.jl")
include("../utilities.jl")
include("gradients_psi_Haux.jl")


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

    # Read data
    #psi = deserialize("../TEMP_psi_v2.dat") # for system_defs_02

    update_from_wavefunc!(Ham, psi) # update the potential
    
    # Prepare Haux
    #Haux = psi' * (Ham*psi) * hx # Hsub, subspace Hamiltonian
    # Using diagonal Haux
    ebands1 = sort(randn(Nstates))
    Haux = diagm(0 => ebands1)
    
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

    α = 2e-3 #3e-5
    α_Haux = 0.1

    Urot2 = zeros(Nstates,Nstates)
    E2 = 0.0
    Udagger = zeros(Nstates,Nstates)

    for i in 1:2000
    
        println("----------------")
        println("Enter step: ", i)
        println("----------------")

        d[:,:] = -Kg[:,:]
        d_Haux[:,:] = -Kg_Haux[:,:]

        psi[:,:] = psi[:,:] + α*d[:,:]
        Haux[:,:] = Haux[:,:] + α_Haux*d_Haux[:,:]
    
        # Orthonormalize
        #ortho_sqrt!(psi)
        
        # Orthonormalize (involves rotation)
        Udagger[:,:] = inv(sqrt(psi'*psi)) ./ sqrt(hx)
        psi[:,:] = psi*Udagger

        println("Check ortho 1:")
        display(psi' * psi * hx)

        # Also rotate Haux
        Haux[:,:] = Udagger' * Haux * Udagger

        println("Haux before (after rotated by Udagger):")
        display(Haux); println()
        
        Urot2[:,:] = transform_psi_Haux!(psi, Haux)

        println("Check ortho 2:")
        display(psi' * psi * hx)
        
        println("Haux after (should be diagonal):")
        display(Haux); println()
    
        E2 = calc_Lfunc_Haux!(Ham, psi, Haux)
        calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)
        
        prec_invK!(Ham, g, Kg)

        println("dot(g,g) = ", dot(g,g)*hx)

        println("g_Haux = ")
        display(g_Haux)

        # why????
        #g_Haux[:,:] = Urot2 * g_Haux * Urot2'
        #Kg_Haux[:,:] = Urot2 * Kg_Haux * Urot2'

        #g[:,:] = g*Urot2
        #Kg[:,:] = Kg*Urot2

        println("E1 = ", E1)
        println("E2 = ", E2)

        if E1 < E2
            println("WARNING: Energy does not decrease")
        end

        dg = dot(g_Haux,g_Haux)*hx
        println("dot g_Haux,g_Haux = ", dg)

        if (abs(E1 - E2) < 1e-7) && (dg < 1e-7)
            println("Converged")
            break
        end

        E1 = E2 # set old value
        Urot[:,:] = Urot2[:,:]

    end

    println("ebands = ")
    display(Ham.electrons.ebands); println()
    
    println("Focc = ")
    println(Ham.electrons.Focc); println()

end
main()