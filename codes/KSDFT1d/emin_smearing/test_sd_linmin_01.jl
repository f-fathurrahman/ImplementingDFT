#  need to run setup_path first

function linmin_quad_v01!(
    α_t,
    Ham, psis, Haux, Hsub, g, g_Haux,
    Kg, Kg_Haux,
    d, d_Haux, rots_cache,
    E_old
)

    gd = 2*dot(g,d) + dot(g_Haux, d_Haux)
    println("gd = $(gd)")
    if gd > 0
        error("Bad step direction")
    end

    NtryMax = 5

    α_safe =  1e-5 # safe step size
    α = α_safe # declare this outside for loop, set to a "safe" value
    α_prev = 0.0
    α_t_ReduceFactor = 0.1
    α_t_IncreaseFactor = 3.0
    is_success = false
    for itry in 1:NtryMax
        #
        println("--- Begin itry linmin trial step = $(itry) using α_t=$(α_t)")
        #
        do_step_psis_Haux!(α_t - α_prev, Ham, psis, Haux, d, d_Haux, rots_cache)
        #
        α_prev = α_t
        update_from_ebands!(Ham)
        update_from_psis!(Ham, psis)
        E_t = calc_Lfunc(Ham, psis)        
        #
        if !isfinite(E_t)
            α_t *= α_t_ReduceFactor
            println("α_t is reduced to=$(α_t)")
            continue # continue
        end
        # prediciton of step size
        c = ( E_t - (E_old + α_t*gd) ) / α_t^2
        α = -gd/(2*c)
        if α < 0
            println("Wrong curvature, α is negative: E_t=$(E_t), E_old=$(E_old)")
            α_t *= α_t_IncreaseFactor
            println("Trial step will become true step. α_t will be set to $(α_t)")
            # calculate gradients
            calc_grad!(Ham, psis, g, Kg, Hsub)
            calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
            # return trial energy and status
            is_success = true
            return E_t, is_success, α_t
        end
        break
    end
    println("Find α = $(α)")
    
    # actual step
    for itry in 1:NtryMax
        #
        println("--- Begin itry linmin actual step = $(itry) using α=$(α)")
        #
        do_step_psis_Haux!(α - α_prev, Ham, psis, Haux, d, d_Haux, rots_cache)
        α_prev = α
        # calculate energy and gradients
        update_from_ebands!(Ham)
        update_from_psis!(Ham, psis)
        E_t2 = calc_Lfunc(Ham, psis)
        calc_grad!(Ham, psis, g, Kg, Hsub)
        calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
        #
        println("Actual step energy 2: E_t2 = $(E_t2)")
        #
        if !isfinite(E_t2)
            α *= α_t_ReduceFactor
            println("α is reduced to=$(α)")
        end
        # trial energy is higher, reduce α
        if E_t2 > E_old
            α *= α_t_ReduceFactor
            println("Energy is not decreasing, try do decrease α to $(α)")
            continue # continue iteration
        else
            println("Actual step is successful")
            is_success = true
            return E_t2, is_success, α
        end
    end

    # default is unsuccessful try
    return Inf, false, α

end



function test_sd_linmin(Ham; NiterMax=100, psis=nothing, Haux=nothing)

    hx = Ham.grid.hx
    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin
    Npoints = Ham.grid.Npoints

    Random.seed!(1234)
    if isnothing(psis)
        psis = Vector{Matrix{Float64}}(undef, Nspin)
        for ispin in 1:Nspin
            psis[ispin] = generate_random_wavefunc(Ham)
        end
    end

    if isnothing(Haux)
        Haux = Vector{Matrix{Float64}}(undef, Nspin)
        for ispin in 1:Nspin
            Haux[ispin] = randn(Float64, Nstates, Nstates)
            Haux[ispin] = 0.5*(Haux[ispin] + Haux[ispin]')
        end
    end

    g = Vector{Matrix{Float64}}(undef, Nspin)
    Kg = Vector{Matrix{Float64}}(undef, Nspin)
    Hsub = Vector{Matrix{Float64}}(undef, Nspin)
    g_Haux = Vector{Matrix{Float64}}(undef, Nspin)
    Kg_Haux = Vector{Matrix{Float64}}(undef, Nspin)
    d = Vector{Matrix{Float64}}(undef, Nspin)
    d_Haux = Vector{Matrix{Float64}}(undef, Nspin)    
    for ispin in 1:Nspin
        g[ispin] = zeros(Float64, Npoints, Nstates)
        Kg[ispin] = zeros(Float64, Npoints, Nstates)
        #
        Hsub[ispin] = zeros(Float64, Nstates, Nstates)
        g_Haux[ispin] = zeros(Float64, Nstates, Nstates)
        Kg_Haux[ispin] = zeros(Float64, Nstates, Nstates)
        #
        g[ispin] = zeros(Float64, Npoints, Nstates)
        g_Haux[ispin] = zeros(Float64, Nstates, Nstates)
    end

    rots_cache = RotationsCache(Nspin, Nstates)

    transform_psis_Haux_update_ebands!(Ham, psis, Haux, rots_cache)
    update_from_ebands!(Ham)
    update_from_psis!(Ham, psis)
    
    E1 = calc_Lfunc(Ham, psis)
    println("Initial energy, E1 = ", E1)

    calc_grad!(Ham, psis, g, Kg, Hsub)
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
    #
    println("test grad psis = ", 2*dot(psis,g)*hx)
    println("test grad Haux = ", dot(Haux,g_Haux))
    # rotate gradients
    rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)

    α_t_start = 1.0
    α_t_min = 1e-6
    α_t = α_t_start

    ebands = Ham.electrons.ebands
    Focc = Ham.electrons.Focc

    for iterSD in 1:NiterMax

        println("\nBegin iterSD = ", iterSD)

        # Set direction
        for ispin in 1:Nspin
            d[ispin] = -Kg[ispin]
            d_Haux[ispin] = -Kg_Haux[ispin]
            constrain_search_dir!(d[ispin], psis[ispin], hx)
        end

        #
        # Do line minimization:
        E_new, is_success, α = linmin_quad_v01!(
            α_t,
            Ham, psis, Haux, Hsub, g, g_Haux, Kg, Kg_Haux, d, d_Haux, rots_cache, E1
        )
        println("Test grad psis before rotate: $(2*dot(g, psis))")
        println("Test grad Haux before rotate: $(dot(Haux, g_Haux))")
        rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)
 
        #
        if is_success
            α_t = α
            println("linminQuad is successful. α_t is updated to α = ", α)
            if α_t < α_t_min
                # bad step size: make sure next test step size is not too bad
                α_t = α_t_start 
                println("Bad step size is encountered, α_t is set to α_t_start = ", α_t_start)
            end
        else
            @warn "Line minimization is not successful"
        end

        @printf("Eigenvalues\n")
        for ist in 1:Nstates
            @printf("%5d %18.10f occ=%10.5f\n", ist, ebands[ist,1], Focc[ist,1])
        end

        dE = abs(E_new - E1)
        @printf("iterSD: %4d %18.10f %10.5e\n", iterSD, E_new, dE)
        if dE < 1e-6
            println("!!!! CONVERGED !!!!!")
            break
        end

        E1 = E_new

    end

    serialize("psis.jldat", psis)
    serialize("Haux.jldat", Haux)

    #@infiltrate

end