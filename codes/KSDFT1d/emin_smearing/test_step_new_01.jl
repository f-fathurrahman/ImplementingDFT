#  need to run setup_path first

function test_step_new_01(Ham)

    hx = Ham.grid.hx
    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin
    Npoints = Ham.grid.Npoints

    Random.seed!(1234)
    psis = Vector{Matrix{Float64}}(undef, Nspin)
    for ispin in 1:Nspin
        psis[ispin] = generate_random_wavefunc(Ham)
    end

    Haux = Vector{Matrix{Float64}}(undef, Nspin)
    for ispin in 1:Nspin
        Haux[ispin] = randn(Float64, Nstates, Nstates)
        Haux[ispin] = 0.5*(Haux[ispin] + Haux[ispin]')
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
    # rotate gradients
    rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)

    # Set direction
    d = -Kg
    d_Haux = -Kg_Haux

    gd = 2*dot(g,d) + dot(g_Haux,d_Haux)
    println("gd = ", gd)
    if gd > 0
        error("Bad step direction")
    end

    # Move forward, positive α
    α = 1.0
    do_step_psis_Haux!(α, Ham, psis, Haux, d, d_Haux, rots_cache)
    update_from_ebands!(Ham)
    update_from_psis!(Ham, psis)
    Etot2 = calc_Lfunc(Ham, psis)
    println("Etot2 = ", Etot2)

    # Move backward, negative α
    α = -1.0
    do_step_psis_Haux!(α, Ham, psis, Haux, d, d_Haux, rots_cache)
    update_from_ebands!(Ham)
    update_from_psis!(Ham, psis)
    Etot3 = calc_Lfunc(Ham, psis)
    println("Etot3 = ", Etot3)
    println("Etot3 should be the same as Etot1 = ", Etot1)

    @infiltrate

end

