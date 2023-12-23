push!(LOAD_PATH, "../")

import Random
using Printf
using LinearAlgebra
using Serialization

using KSDFT1d

#include("system_defs_01.jl")
include("system_defs_02.jl")
#include("system_defs_03.jl")

include("Lfunc.jl")
include("../utilities.jl")
include("gradients_psi_Haux.jl")

include("linemin_grad.jl")
include("linemin_armijo.jl")
include("linemin_quad.jl")


function solve_Emin_SD!(Ham, psi, Haux, g, g_Haux, Kg, Kg_Haux, d, d_Haux)

    hx = Ham.grid.hx
    Npoints = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates

    Hsub = zeros(Float64, size(Haux))

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
    α = 1e-4
    α_Haux = 1e-1
    β = 0.0
    β_Haux = 0.0

    is_linmin_success1 = true
    is_linmin_success2 = true

    g_old = zeros(Float64, size(psi))
    Kg_old = zeros(Float64, size(psi))
    d_old = zeros(Float64, size(psi))

    g_Haux_old = zeros(Float64, size(Haux))
    Kg_Haux_old = zeros(Float64, size(Haux))
    d_Haux_old = zeros(Float64, size(Haux))


    println("E initial = ", E1)

    for iterEmin in 1:500

        println("\nBegin iterEmin = ", iterEmin)

#=
        if iterEmin >= 2
            num1 = 2*dot(g-g_old, Kg)*hx
            denum1 = 2*dot(g_old, Kg_old)*hx
            β = num1/denum1
            if β < 0.0
                β = 0.0
            end
            num2 = dot(g_Haux-g_Haux_old, Kg_Haux)
            denum2 = dot(g_Haux_old, Kg_Haux_old)
            β_Haux = num2/denum2
            if β_Haux < 0.0
                β_Haux = 0.0
            end
        end
=#

        println("β = $(β) , β_Haux = $(β_Haux)")

        # Search direction
        d[:,:] = -Kg + β*d_old
        d_Haux[:,:] = -Kg_Haux + β_Haux*d_Haux_old
        #
        constrain_search_dir!(d, psi, hx)

        α, is_linmin_success1 = linemin_armijo_psi(Ham, psi, Haux, d, d_Haux, E1, α0=1.0, reduce_factor=0.25)
        α_Haux, is_linmin_success2 = linemin_armijo_Haux(Ham, psi, Haux, d, d_Haux, E1, α0=1.0, reduce_factor=0.25)
        
        #α, is_linmin_success1 = linemin_quad_psi(Ham, psi, Haux, g, g_Haux, d, d_Haux, E1)
        #α_Haux, is_linmin_success2 = linemin_quad_Haux(Ham, psi, Haux, g, g_Haux, d, d_Haux, E1)

        println("α = ", α, " α_Haux = ", α_Haux)

        # We will stop iteration if line minimization is not successful
        if !(is_linmin_success1 && is_linmin_success2)
            is_converged = false
            break
        end

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
        # Update psi
        psi[:,:] = psi + α*d
        Haux[:,:] = Haux + α_Haux*d_Haux
        prepare_psi_Haux!(psi, Haux, hx)
        # Evaluate
        E1 = calc_Lfunc_Haux!(Ham, psi, Haux)
        calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)
        #
        calc_grad_no_Focc!(Ham, psi, Kg)
        prec_invK!(Ham, Kg)
        #
        dE = E1 - E_old
        dg_Haux = dot(g_Haux,g_Haux)*hx/length(g_Haux)
        dg = dot(g,g)*hx/length(g)
        #
        @printf("EMIN: %8d %18.10f %18.10e [%18.10e,%18.10e]\n", iterEmin, E1, dE, dg, dg_Haux)
        if E1 > E_old
            println("!!! WARNING: Energy is not decreasing")
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
    #iseed = abs(rand(Int64))
    iseed = 1234
    println("iseed = ", iseed)
    Random.seed!(iseed)
    #Random.seed!(1) # vary this to find problematic case?

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
        prepare_psi_Haux!(psi, Haux, hx)
    end

    println("ebands = ")
    display(Ham.electrons.ebands); println()
    
    println("Focc = ")
    println(Ham.electrons.Focc); println()

#    return
end

main()
