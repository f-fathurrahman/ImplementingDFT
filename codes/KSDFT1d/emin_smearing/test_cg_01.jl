#  need to run setup_path first

function test_cg_01(Ham; NiterMax=100, psis=nothing, Haux=nothing)

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
    gPrev = Vector{Matrix{Float64}}(undef, Nspin)
    gPrev_Haux = Vector{Matrix{Float64}}(undef, Nspin)  
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
        #
        d[ispin] = zeros(Float64, Npoints, Nstates)
        d_Haux[ispin] = zeros(Float64, Nstates, Nstates)
        #
        gPrev[ispin] = zeros(Float64, Npoints, Nstates)
        gPrev_Haux[ispin] = zeros(Float64, Nstates, Nstates)
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

    do_force_grad_dir = true
    gKNorm = 0.0
    gKNormPrev = 0.0 # current and previous norms of the preconditioned gradient

    for iterCG in 1:NiterMax

        println("\nBegin iterCG = ", iterCG)

        gKNorm = 2*dot(g, Kg)*hx + dot(g_Haux, Kg_Haux)

        β = 0.0
        if !do_force_grad_dir
            gd = 2*dot(g, d)*hx + dot(g_Haux, d_Haux)
            gPrevKg = 2*dot(gPrev, Kg)*hx + dot(gPrev_Haux, Kg_Haux)
            gg = 2*dot(g, g)*hx + dot(g_Haux, g_Haux)
            dd = 2*dot(d, d)*hx + dot(d_Haux, d_Haux)
            @printf("linmin: %10.3le", gd/sqrt(gg*dd))
            @printf("  cgtest: %10.3le\n", gPrevKg/sqrt(gKNorm*gKNormPrev))
            # Update beta:
            β = (gKNorm - gPrevKg)/gKNormPrev
            if β < 0.0
                println("Resetting CG")
            end
            println("β = ", β)                    
        end

        do_force_grad_dir = false

        # XXX TODO Check convergence here?

        # Save previous gradient
        gKNormPrev = gKNorm
        for ispin in 1:Nspin
            gPrev[ispin] = g[ispin]
            gPrev_Haux[ispin] = g_Haux[ispin]
        end

        # Set direction
        for ispin in 1:Nspin
            d[ispin] = -Kg[ispin] + β*d[ispin]
            d_Haux[ispin] = -Kg_Haux[ispin] + β*d_Haux[ispin]
            constrain_search_dir!(d[ispin], psis[ispin], hx)
        end    

        #
        # Do line minimization:
        E_new, is_success, α = linmin_quad_v01!(
            α_t,
            Ham, psis, Haux, Hsub, g, g_Haux, Kg, Kg_Haux, d, d_Haux, rots_cache, E1
        )
        println("test grad psis = ", 2*dot(psis,g)*hx)
        println("test grad Haux = ", dot(Haux,g_Haux))
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
            # XXX Reset search direction here ?
            # do_force_grad_dir set to true ?
        end

        @printf("Eigenvalues\n")
        for ist in 1:Nstates
            @printf("%5d %18.10f occ=%10.5f\n", ist, ebands[ist,1], Focc[ist,1])
        end

        dE = abs(E_new - E1)
        @printf("iterCG: %4d %18.10f %10.5e\n", iterCG, E_new, dE)
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