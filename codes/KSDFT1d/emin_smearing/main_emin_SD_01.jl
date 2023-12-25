push!(LOAD_PATH, "../")

import Random
using Printf
using LinearAlgebra
using Serialization

using KSDFT1d

include("system_defs_01.jl")
#include("system_defs_02.jl")
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

    E_old = Inf
    dE = Inf
    dg_Haux = Inf
    dg = Inf

    Nconverges = 0
    is_increasing = false
    is_converged = false
    α = 1.0
    β = 0.0

    g_old = zeros(Float64, size(psi))
    Kg_old = zeros(Float64, size(psi))
    d_old = zeros(Float64, size(psi))

    g_Haux_old = zeros(Float64, size(Haux))
    Kg_Haux_old = zeros(Float64, size(Haux))
    d_Haux_old = zeros(Float64, size(Haux))

    # Evaluate total energy and gradient by calling Lfunc
    E1 = calc_Lfunc_Haux!(Ham, psi, Haux)
    calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)    
    calc_grad_no_Focc!(Ham, psi, Kg)
    prec_invK!(Ham, Kg)

    println("Initial E = ", E1)


    for iterEmin in 1:500

        println("\nBegin iterEmin = ", iterEmin)


        if iterEmin >= 2
            #num1 = 2*dot(g-g_old, Kg)*hx + dot(g_Haux-g_Haux_old, Kg_Haux)
            #denum1 = 2*dot(g_old, Kg_old)*hx + dot(g_Haux_old, Kg_Haux_old)
            num1 = 2*dot(g,Kg)*hx + dot(g_Haux,Kg_Haux)
            denum1 = 2*dot(g_old,Kg_old)*hx + dot(g_Haux_old,Kg_Haux_old)
            β = num1/denum1
            if β < 0.0
                β = 0.0
            end
        end
        #println("β = $(β) , β_Haux = $(β_Haux)")

        println("β = $(β)")


        # Save old search direction
        d_old[:] .= d[:]
        d_Haux_old[:] .= d_Haux[:]
        #
        # Set search direction
        #
        d[:,:] = -Kg + β*d_old
        d_Haux[:,:] = -Kg_Haux + β*d_Haux_old
        constrain_search_dir!(d, psi, hx)

        gd_psi = 2*dot(g, d)*hx
        gd_Haux = dot(g_Haux, d_Haux)
        gdotd = gd_psi + gd_Haux
        println("gd_psi = ", gd_psi)
        println("gd_Haux = ", gd_Haux)
        println("gdotd = ", gdotd)
        #
        if gdotd > 0
            println("CG: !!! Bad step direction, reset CG")
            d[:,:] = -Kg
            d_Haux[:,:] = -Kg_Haux
            constrain_search_dir!(d, psi, hx)
            gd_psi = 2*dot(g, d)*hx
            gd_Haux = dot(g_Haux, d_Haux)
            gdotd = gd_psi + gd_Haux
            println("gd_psi = ", gd_psi)
            println("gd_Haux = ", gd_Haux)
            println("gdotd = ", gdotd)
        end


        # Line minimization
        α, is_linmin_success = linemin_armijo(Ham, psi, Haux, d, d_Haux, E1)

        #α = linemin_grad(Ham, psi, Haux, g, g_Haux, d, d_Haux, E1)
        #is_linmin_success = true # force to true

        #α, is_linmin_success = linemin_quad(Ham, psi, Haux, g, g_Haux, d, d_Haux, E1)
        #println("α = ", α)


        num1 = 2*dot(g,d_old)*hx + dot(g_Haux,d_Haux_old)
        gg = 2*dot(g,g)*hx + dot(g_Haux,g_Haux)
        dd = 2*dot(d_old,d_old)*hx + dot(d_Haux_old,d_Haux_old)
        println("num1 = ", num1)
        println("gg*dd = ", gg*dd)
        linmin_test = num1/sqrt(gg*dd)
        println("linmin_test = ", linmin_test)


        #
        # Save old variables
        #
        E_old = E1
        g_old[:] .= g[:]
        Kg_old[:] .= Kg[:]
        g_Haux_old[:] .= g_Haux[:]
        Kg_Haux_old[:] .= Kg_Haux[:]
        #
        # Actual step
        #
        psi[:,:] = psi + α*d
        Haux[:,:] = Haux + α*d_Haux
        prepare_psi_Haux!(psi, Haux, hx)
        #
        E1 = calc_Lfunc_Haux!(Ham, psi, Haux)
        calc_grad_Lfunc_Haux!(Ham, psi, Haux, g, Hsub, g_Haux, Kg_Haux)
        calc_grad_no_Focc!(Ham, psi, Kg)
        prec_invK!(Ham, Kg)
        #
        #
        dE = E1 - E_old
        dg_Haux = dot(g_Haux,g_Haux)/length(g_Haux)
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
        # Back to last search ?
        psi[:,:] = psi[:,:] + g
        Haux[:,:] = Haux[:,:] + g_Haux
        prepare_psi_Haux!(psi, Haux, hx)
    end

    println("ebands = ")
    display(Ham.electrons.ebands); println()
    
    println("Focc = ")
    println(Ham.electrons.Focc); println()

#    return
end

main()
