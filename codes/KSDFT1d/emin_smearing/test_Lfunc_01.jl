#  need to run setup_path first

function test_Lfunc_01(Ham)

    hx = Ham.grid.hx
    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin

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

    rots_cache = RotationsCache(Nspin, Nstates)
    transform_psis_Haux_update_ebands!(Ham, psis, Haux, rots_cache)

    update_from_ebands!(Ham)
    update_from_psis!(Ham, psis)
    Etot1 = calc_Lfunc(Ham, psis)

    println("Etot1 = ", Etot1)

    @infiltrate

end

#=
# Check invariance w.r.t unitary transform
println("\nUsing rotated psi")
psi2 = psi*Urot
Etot2 = calc_Lfunc_ebands!(Ham, psi*Urot, ebands)
=#