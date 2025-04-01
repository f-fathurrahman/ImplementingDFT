#  need to run setup_path first

function test_grad_new_01(Ham)

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
    for ispin in 1:Nspin
        g[ispin] = zeros(Float64, Npoints, Nstates)
        Kg[ispin] = zeros(Float64, Npoints, Nstates)
        Hsub[ispin] = zeros(Float64, Nstates, Nstates)
        g_Haux[ispin] = zeros(Float64, Nstates, Nstates)
        Kg_Haux[ispin] = zeros(Float64, Nstates, Nstates)
    end

    rots_cache = RotationsCache(Nspin, Nstates)

    psis_orig = copy(psis)
    Haux_orig = copy(Haux)

    transform_psi_Haux_update_ebands!(Ham, psis, Haux, rots_cache)
    update_from_ebands!(Ham)
    update_from_psis!(Ham, psis)
    
    Etot1 = calc_Lfunc(Ham, psis)
    println("Etot1 = ", Etot1)

    calc_grad!(Ham, psis, g, Kg, Hsub)
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
    # rotate gradients
    rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)

    @infiltrate

end
