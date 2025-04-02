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