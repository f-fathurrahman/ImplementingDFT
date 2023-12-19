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

function linemin_armijo(Ham, psi, Haux, d_in, d_Haux_in, E1;
    α0=10.0, reduce_factor=0.5, α_safe=1e-8
)
    #
    hx = Ham.grid.hx
    #
    τ = reduce_factor # reduction factor
    #
    α = α0
    d = α*d_in
    d_Haux = α*d_Haux_in
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
    for iterE in 1:20
        #
        dE = E_new - E1
        @printf("LineMin: %3d %18.5e %18.5e\n", iterE, α, dE)
        if dE < 0.0
            println("E_new is smaller, α=$(α) is accepted, found in iterE=$(iterE)")
            is_success = true
            break
        end        
        #
        α = τ*α
        d[:,:] = α*d_in
        d_Haux[:,:] = α*d_Haux_in
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
        return α, true
    else
        #println("LineMin: unsuccessful, returning last α")
        #return α
        #error("LineMin: unsuccessful")
        println("LineMin: unsuccessful, returning α_safe")
        return α_safe, true #false
    end
end



function solve_Emin_SD!(Ham, psi, Haux, g, g_Haux, Kg, Kg_Haux, d, d_Haux)

    hx = Ham.grid.hx
    Npoints = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates

    #
    Hsub = zeros(Float64, size(Haux))
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
    is_converged = false
    α = 1.0
    β = 0.0
    β_Haux = 0.0


    g_old = zeros(Float64, size(psi))
    Kg_old = zeros(Float64, size(psi))
    d_old = zeros(Float64, size(psi))

    g_Haux_old = zeros(Float64, size(Haux))
    Kg_Haux_old = zeros(Float64, size(Haux))
    d_Haux_old = zeros(Float64, size(Haux))


    println("E initial = ", E1)

    for iterEmin in 1:2000


        if iterEmin >= 3
            #β = real( dot(g - g_old, Kg) )/real( dot(g_old,Kg_old) )
            numb = dot(g - g_old, Kg) + dot(g_Haux - g_Haux_old, Kg_Haux)
            denumb = dot(g_old,Kg_old) + dot(g_Haux_old,Kg_Haux_old)
            β = numb/denumb
            if β < 0.0
                β = 0.0
            end
            #β_Haux = real( dot(g_Haux - g_Haux_old, Kg_Haux) )/real( dot(g_Haux_old,Kg_Haux_old) )
            #if β_Haux < 0.0
            #    β_Haux = 0.0
            #end
        end
        #println("β = $(β) , β_Haux = $(β_Haux)")
        println("β = $(β)")


        #d[:,:] = -Kg + (β+β_Haux)*d_old
        #d_Haux[:,:] = -Kg_Haux + (β+β_Haux)*d_Haux_old

        d[:,:] = -Kg + β*d_old
        d_Haux[:,:] = -Kg_Haux + β*d_Haux_old


        println()
        if is_increasing
            α, is_linmin_success = linemin_armijo(Ham, psi, Haux, d, d_Haux, E1, reduce_factor=0.1)
        else
            α, is_linmin_success = linemin_armijo(Ham, psi, Haux, d, d_Haux, E1)
        end

        if !is_linmin_success
            is_converged = false
            break
        end

        #
        #println("DEBUG: dot grad psi: ", dot(g,Kg)*hx)
        #println("DEBUG: dot grad Haux: ", dot(g_Haux,Kg_Haux))
        #
        #println("DEBUG: dot step psi: ", dot(g,d)*hx)
        #println("DEBUG: dot step Haux: ", dot(g_Haux,d_Haux))
        #
        # Save old variables
        E_old = E1
        #
        g_old[:] .= g[:]
        Kg_old[:] .= Kg[:]
        d_old[:] .= d[:]
        #
        g_Haux_old[:] .= g_Haux[:]
        Kg_Haux_old[:] .= Kg_Haux[:]
        d_Haux_old[:] .= d_Haux[:]
        #
        #
        psi[:,:] = psi + α*d
        Haux[:,:] = Haux + α*d_Haux
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
            is_increasing = true
            is_converged = false
            break
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
            is_converged = true
            break
        end
        #
    end
    return E1, is_converged, α
end


function main()
    Ham = init_Hamiltonian()

    hx = Ham.grid.hx
    Npoints = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates

    # Random wavefunc
    #Random.seed!(111) # vary this to find problematic case?
    psi = generate_random_wavefunc(Ham)
    Haux = diagm( 0 => sort(randn(Nstates)) )
    #
    d = zeros(Float64, size(psi))
    d_Haux = zeros(Float64, size(Haux))
    # Allocate memory
    g = zeros(Float64, size(psi))
    Kg = zeros(Float64, size(psi))
    #
    g_Haux = zeros(Float64, size(Haux))
    Kg_Haux = zeros(Float64, size(Haux))


    Etot = Inf
    is_converged = false
    α = 1.0

    NTRY_MAX = 1
    for itry in 1:NTRY_MAX
        Etot, is_converged, α = solve_Emin_SD!(Ham, psi, Haux, g, g_Haux, Kg, Kg_Haux, d, d_Haux)
        if is_converged
            println("Success in itry=$(itry)")
            break
        end
        #
        psi[:,:] = psi[:,:] + g
        Haux[:,:] = Haux[:,:] + g_Haux
        #
        Udagger = inv(sqrt(psi'*psi)) ./ sqrt(hx)
        psi[:,:] = psi*Udagger
        Haux[:,:] = Udagger' * Haux * Udagger
        Urot = transform_psi_Haux!(psi, Haux)
    end

    println("ebands = ")
    display(Ham.electrons.ebands); println()
    
    println("Focc = ")
    println(Ham.electrons.Focc); println()

#    return
end

main()
