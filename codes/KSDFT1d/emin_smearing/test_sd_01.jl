#  need to run setup_path first

function test_sd_01(Ham; NiterMax=100, psis=nothing, Haux=nothing)

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

    psis_orig = copy(psis)
    Haux_orig = copy(Haux)

    transform_psis_Haux_update_ebands!(Ham, psis, Haux, rots_cache)
    update_from_ebands!(Ham)
    update_from_psis!(Ham, psis)
    
    Etot1 = calc_Lfunc(Ham, psis)
    println("Etot1 = ", Etot1)

    calc_grad!(Ham, psis, g, Kg, Hsub)
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
    #
    println("test grad psis = ", 2*dot(psis,g)*hx)
    println("test grad Haux = ", dot(Haux,g_Haux))
    # rotate gradients
    rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)

    α = 0.01
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

        gd = 2*dot(g,d)*hx + dot(g_Haux,d_Haux)
        gg = 2*dot(g,g)*hx + dot(g_Haux,g_Haux)
        println("gd = ", gd, " gg = ", gg)
        if gd > 0
            error("Bad step direction")
        end

        do_step_psis_Haux!(α, Ham, psis, Haux, d, d_Haux, rots_cache)
        update_from_ebands!(Ham)
        update_from_psis!(Ham, psis)
        Etot_new = calc_Lfunc(Ham, psis)
        #
        calc_grad!(Ham, psis, g, Kg, Hsub)
        calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
        println("test grad psis = ", 2*dot(psis,g)*hx)
        println("test grad Haux = ", dot(Haux,g_Haux))
        # rotate gradients
        rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)
        #
        if Etot_new > Etot1
            #
            error("Higher Etot_new")
            #
            # undo step
            do_step_psis_Haux!(-α, Ham, psis, Haux, d, d_Haux, rots_cache)
            update_from_ebands!(Ham)
            update_from_psis!(Ham, psis)
            #
            Etot1 = calc_Lfunc(Ham, psis)
            #
            calc_grad!(Ham, psis, g, Kg, Hsub)
            calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
            # rotate gradients
            rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)
            # Reduce α
            α = 0.5*α
            println("WARN: Etot_new is higher, reducing α to ", α)
            println("Should return to prev points: Etot1 = ", Etot1)
            continue # next iter
        end

        @printf("Eigenvalues\n")
        for ist in 1:Nstates
            @printf("%5d %18.10f occ=%10.5f\n", ist, ebands[ist,1], Focc[ist,1])
        end

        dE = abs(Etot_new - Etot1)
        @printf("%4d %18.10f %10.5e\n", iterSD, Etot_new, dE)
        if dE < 1e-6
            println("!!!! CONVERGED !!!!!")
            break
        end

        Etot1 = Etot_new

    end

    serialize("psis.jldat", psis)
    serialize("Haux.jldat", Haux)

    @infiltrate

end

