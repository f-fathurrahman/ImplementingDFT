push!(LOAD_PATH, "../")

import Random
using Printf
using LinearAlgebra
using Serialization

using KSDFT1d

#include("system_defs_01.jl")
#include("system_defs_02.jl")
include("system_defs_03.jl")

include("Lfunc.jl")
include("../utilities.jl")
include("gradients_psi_Haux.jl")

function linemin_armijo(Ham, psi, Haux, Kg, Kg_Haux, E1;
    α0=10.0, reduce_factor=0.5, α_safe=1e-8
)
    #
    hx = Ham.grid.hx
    #
    τ = reduce_factor # reduction factor
    #
    α = α0
    d = -α*Kg
    d_Haux = -α*Kg_Haux
    #
    psi_new = psi + d
    Haux_new = Haux + d_Haux
    #
    Udagger = inv(sqrt(psi_new'*psi_new)) ./ sqrt(hx)
    psi_new[:,:] = psi_new*Udagger
    Haux_new = Udagger' * Haux_new * Udagger
    Urot = transform_psi_Haux!(psi_new, Haux_new)
    #
    E_new = calc_Lfunc_Haux!(Ham, psi_new, Haux_new)
    #
    is_success = false
    for iterE in 1:40
        #
        dE = E_new - E1
        @printf("LineMin: %3d %18.5e %18.5e\n", iterE, α, dE)
        if dE < 0.0
            println("E_new is smaller, α=$(α) is accepted")
            is_success = true
            break
        end        
        #
        α = τ*α
        d[:,:] = -α*Kg
        d_Haux[:,:] = -α*Kg_Haux
        #
        psi_new[:,:] = psi + d
        Haux_new[:,:] = Haux + d_Haux
        #
        Udagger[:,:] = inv(sqrt(psi_new'*psi_new)) ./ sqrt(hx)
        psi_new[:,:] = psi_new*Udagger
        Haux_new[:,:] = Udagger' * Haux_new * Udagger
        Urot[:,:] = transform_psi_Haux!(psi_new, Haux_new)
        #
        E_new = calc_Lfunc_Haux!(Ham, psi_new, Haux_new)
    end

    if is_success
        return α
    else
        #println("LineMin: unsuccessful, returning last α")
        #return α
        #error("LineMin: unsuccessful")
        println("LineMin: unsuccessful, returning α_safe")
        return α_safe
    end
end



#function main()


    Ham = init_Hamiltonian()

    hx = Ham.grid.hx
    Npoints = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates

    # Random wavefunc
    Random.seed!(1) # vary this to find problematic case?
    psi = generate_random_wavefunc(Ham)
    Haux = diagm( 0 => sort(randn(Nstates)) )

    # psi and Haux are ready

    # Allocate memory
    g = zeros(Float64, size(psi))
    Kg = zeros(Float64, size(psi))
    #
    g_Haux = zeros(Float64, size(Haux))
    Kg_Haux = zeros(Float64, size(Haux))
    #
    Hsub = zeros(Float64, size(Haux))
    #
    d = zeros(Float64, size(psi))
    d_Haux = zeros(Float64, size(Haux))
    #
    Udagger = zeros(Float64, size(Haux))
    Urot = zeros(Float64, size(Haux))

    # Evaluate total energy and gradient by calling Lfunc
    E1 = calc_Lfunc_Haux!(Ham, psi, Haux)
    calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)
    
    #prec_invK!(Ham, g, Kg) # Precondition
    calc_grad_no_Focc!(Ham, psi, Kg)
    prec_invK!(Ham, Kg)


    E_old = E1

    dE = Inf
    dg_Haux = Inf
    dg = Inf

    Nconverges = 0
    is_increasing = false

    for iterEmin in 1:2000

        println()
        if is_increasing
            α = linemin_armijo(Ham, psi, Haux, Kg, Kg_Haux, E1, reduce_factor=0.1)
        else
            α = linemin_armijo(Ham, psi, Haux, Kg, Kg_Haux, E1)
        end
        d[:,:] = -α*Kg
        d_Haux[:,:] = -α*Kg_Haux
        #
        println("DEBUG: dot grad psi: ", dot(g,Kg)*hx)
        println("DEBUG: dot grad Haux: ", dot(g_Haux,Kg_Haux))
        #
        println("DEBUG: dot step psi: ", dot(g,d)*hx)
        println("DEBUG: dot step Haux: ", dot(g_Haux,d_Haux))
        #
        #
        psi[:,:] = psi + d
        Haux[:,:] = Haux + d_Haux
        #
        Udagger[:,:] = inv(sqrt(psi'*psi)) ./ sqrt(hx)
        psi[:,:] = psi*Udagger
        Haux[:,:] = Udagger' * Haux * Udagger
        Urot[:,:] = transform_psi_Haux!(psi, Haux)
        #
        E1 = calc_Lfunc_Haux!(Ham, psi, Haux)
        calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)
        #prec_invK!(Ham, g, Kg) # Precondition
        calc_grad_no_Focc!(Ham, psi, Kg)
        prec_invK!(Ham, Kg)
        #
        #
        dE = E1 - E_old
        dg_Haux = dot(g_Haux,g_Haux)*hx/length(g_Haux)
        dg = dot(g,g)*hx/length(g)
        #
        @printf("EMIN: %8d %18.10f %18.10e [%18.10e,%18.10e]\n", iterEmin, E1, dE, dg, dg_Haux)
        if E1 > E_old
            println("!!! Energy is not decreasing")
            #error("!!! Energy is not decreasing")
            is_increasing = true
        else
            is_increasing = false
        end
        if (abs(dE) < 1e-6) && (abs(dg_Haux) < 1e-3) && !is_increasing
            Nconverges += 1
        else
            Nconverges = 0
        end
        if Nconverges >= 2
            println("\n*** Converged ***\n")
            break
        end
        #
        E_old = E1
    end

    println("ebands = ")
    display(Ham.electrons.ebands); println()
    
    println("Focc = ")
    println(Ham.electrons.Focc); println()

#end

#main()
