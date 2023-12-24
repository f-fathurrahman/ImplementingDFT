# Vary psi and Haux simulataneously
function linemin_quad(α_t::Float64, Ham, psi, Haux, g, g_Haux, d, d_Haux, E1)

    hx = Ham.grid.hx

    α_start = 1.0
    α_min = 1e-10
    reduceFactor = 0.1
    increaseFactor = 3.0
    NalphaAdjustMax = 3

    gdotd = 2*dot(g, d)*hx + dot(g_Haux, d_Haux)
    if gdotd > 0.0
        println("linemin_quad: !!! Bad step direction")
        return 1.0, α_t, false
    end

    α_prev = 0.0
    α = 1.0

    for iadjust in 1:NalphaAdjustMax

        if α_t < α_min
            println("linemin_quad: α_t is too small")
            α = α_prev
            return α, α_t, false
        end

        α_prev = α_t
        println("Using α_t = ", α_t)
        psic = psi + α_t*d # trial wavefunc
        Hauxc = Haux + α_t*d_Haux
        #
        prepare_psi_Haux!(psic, Hauxc, hx)
        #
        Etrial = calc_Lfunc_Haux!(Ham, psic, Hauxc)
        ΔEdir = 2*dot(g, α_t*d)*hx + dot(g_Haux, α_t*d_Haux)
        #
        println("linemin_quad: E1     = ", E1)
        println("linemin_quad: Etrial = ", Etrial)
        println("linemin_quad: ΔEdir  = ", ΔEdir)
        println("linemin_quad: ratio of energy diff = ", (Etrial - E1)/ΔEdir)
        #
        c = ( Etrial - (E1 + ΔEdir) ) / α_t
        α = -ΔEdir/(2*c)
        println("linemin_quad: c = ", c)
        println("linemin_quad: α = ", α)
        #if α < 0.0
        if Etrial < E1
            println("linemin_quad: Etrial is lower than E1")
            println("linemin_quad: Wrong curvature, returning α_t")
            #α_t = α_t*increaseFactor
            #println("linemin_quad: α_t is increased to ", α_t)
            return α_t, α_t, true
        end
        #
        if α/α_t > increaseFactor
            α_t *= increaseFactor
            @printf("Predicted α/α_t > %lf, increasing α_t to %le\n", increaseFactor, α_t)
        end
        #
        if α_t/α < reduceFactor
            α_t *= reduceFactor
            @printf("Predicted α/α_t < %lf, reducing α_t to %le\n", reduceFactor, α_t)
        end
        #
        println("linemin_quad is successful: α = ", α, " α_t = ", α_t)
        break
    end

    return α, α_t, true
end