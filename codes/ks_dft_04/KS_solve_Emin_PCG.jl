function calc_grad!(
    Ham::Hamiltonian,
    ψ::Array{Float64,2},
    g::Array{Float64,2},
    Hsub::Matrix{Float64}
)
    ispin = Ham.ispin
    Nstates = size(ψ,2)
    Nspin = Ham.electrons.Nspin

    @views Focc = Ham.electrons.Focc[:,ispin]
    Hψ = op_H( Ham, ψ )
    Hsub[:] = ψ' * Hψ * Ham.grid.dVol
    Hψ = Hψ - ψ*Hsub
    for ist in 1:Nstates
        @views g[:,ist] = Focc[ist] * Hψ[:,ist]
    end

    if Nspin == 1
        need_correction = !all(Focc .== 2)
    else
        need_correction = !all(Focc .== 1)
    end

    println("need_correction = ", need_correction)
    if !need_correction
        # immediate return
        return
    end

    F = Diagonal(Focc)
    HFH = Hsub*F - F*Hsub
    println("HFH = "); display(HFH); println()
    if Nspin == 1
        g[:,:] = g[:,:] + 0.5*ψ*HFH
    else
        g[:,:] = g[:,:] + ψ*HFH
    end
    return
end


# No Hsub
function calc_grad!(
    Ham::Hamiltonian,
    ψ::Array{Float64,2},
    g::Array{Float64,2}
)
    ispin = Ham.ispin
    Nstates = size(ψ,2)
    Nspin = Ham.electrons.Nspin

    @views Focc = Ham.electrons.Focc[:,ispin]
    Hψ = op_H( Ham, ψ )
    Hsub = ψ' * Hψ * Ham.grid.dVol
    Hψ = Hψ - ψ*Hsub
    for ist in 1:Nstates
        @views g[:,ist] = Focc[ist] * Hψ[:,ist]
    end

    if Nspin == 1
        need_correction = !all(Focc .== 2)
    else
        need_correction = !all(Focc .== 1)
    end

    println("need_correction = ", need_correction)
    if !need_correction
        # immediate return
        return
    end

    F = Diagonal(Focc)
    HFH = Hsub*F - F*Hsub
    println("HFH = "); display(HFH); println()
    if Ham.electrons.Nspin == 1
        g[:,:] = g[:,:] + 0.5*ψ*HFH
    else
        g[:,:] = g[:,:] + ψ*HFH
    end
    return
end

function calc_energies_grad!(
    Ham::Hamiltonian,
    psis::Vector{Array{Float64,2}},
    Rhoe::Array{Float64,2},
    g::Vector{Array{Float64,2}},
    Kg::Vector{Array{Float64,2}},
    Hsub::Vector{Matrix{Float64}}
)
    calc_rhoe!( Ham, psis, Rhoe )
    update!( Ham, Rhoe )
    calc_energies!( Ham, psis )

    Nstates = size(psis[1],2)
    Nspin = Ham.electrons.Nspin

    for i in 1:Nspin
        Ham.ispin = i
        calc_grad!( Ham, psis[i], g[i], Hsub[i] )
        # apply preconditioner
        Kg[i][:,:] = g[i][:,:]
        for ist in 1:Nstates
            @views ldiv!(Ham.precKin, Kg[i][:,ist])
        end
    end
    return sum( Ham.energies )
end


function calc_energies_only!( Ham, psis, Rhoe )
    calc_rhoe!( Ham, psis, Rhoe )
    update!( Ham, Rhoe )
    calc_energies!( Ham, psis )
    return sum( Ham.energies )
end


function linmin_grad!(
    Ham::Hamiltonian,
    psis::Vector{Array{Float64,2}},
    g::Vector{Array{Float64,2}},
    d::Vector{Array{Float64,2}},
    Rhoe::Array{Float64,2},
    psisc::Vector{Array{Float64,2}},
    gt::Vector{Array{Float64,2}};
    αt=3e-5
)

    Nspin = Ham.electrons.Nspin
    Nbasis = size(psis[1], 1)
    Nstates = size(psis[1], 2)
    
    dVol= Ham.grid.dVol

    for i in 1:Nspin
        psisc[i][:] = psis[i] + αt*d[i]
        ortho_sqrt!(psisc[i])
        psisc[i][:] = psisc[i][:]/sqrt(dVol)
    end

    calc_rhoe!( Ham, psisc, Rhoe )
    update!( Ham, Rhoe )
    for i in 1:Nspin
        Ham.ispin = i
        calc_grad!( Ham, psisc[i], gt[i] )
    end

    denum = dot(g - gt, d)*dVol

    if denum != 0.0
        α = abs( αt * dot(g, d)*dVol/denum )
    else
        α = 0.0
    end
    println("α in linmin_grad = ", α)
    return α
end

function constrain_search_dir!( d, psi, dVol )
    d[:,:] = d[:,:] - psi * ( psi' * d * dVol )
    return
end

function KS_solve_Emin_PCG!(
    Ham::Hamiltonian, psis::Vector{Array{Float64,2}};
    α_t=3e-5, NiterMax=200,
    etot_conv_thr=1e-6
)
    Nspin = size(psis, 1)
    Nbasis = size(psis[1], 1)
    Nstates = size(psis[1], 2)

    dVol = Ham.grid.dVol

    g = Vector{Array{Float64,2}}(undef,Nspin)
    for i in 1:Nspin
        g[i] = zeros(Float64,Nbasis,Nstates)
    end

    Kg = deepcopy(g)
    gPrev = deepcopy(g)
    d = deepcopy(g)
    psisc = deepcopy(g)
    gt = deepcopy(g)

    Hsub = Vector{Array{Float64,2}}(undef,Nspin)
    for i in 1:Nspin
        Hsub[i] = zeros(Float64,Nstates,Nstates)
    end

    Rhoe = zeros(Float64,Nbasis,Nspin)

    Ham.energies.NN = calc_E_NN( Ham.atoms, Ham.pspots )

    Etot = calc_energies_grad!( Ham, psis, Rhoe, g, Kg, Hsub )

    for i in 1:Nspin
        d[i][:,:] = -Kg[i][:,:]
    end

    #constrain_search_dir!( d, psi, dVol )

    gPrevUsed = true

    α = 0.0
    β = 0.0
    gKnorm = 0.0
    gKnormPrev = 0.0
    force_grad_dir = true

    Etot_old = Etot
    Nconverges = 0
    
    cg_test = 0.0

    @printf("\n")
    @printf("Minimizing Kohn-Sham energy using PCG\n")
    @printf("-------------------------------------\n")
    @printf("NiterMax  = %d\n", NiterMax)
    @printf("α_t       = %e\n", α_t)
    @printf("conv_thr  = %e\n", etot_conv_thr)

    for iter in 1:NiterMax

        gKnorm = dot(g,Kg)*dVol
        if !force_grad_dir
            dotgd = dot(g,d)*dVol
            if gPrevUsed
                dotgPrevKg = dot(gPrev, Kg)*dVol
            else
                dotgPrevKg = 0.0
            end
            β = (gKnorm - dotgPrevKg)/gKnormPrev # Polak-Ribiere
        end

        if β < 0.0
            #println("Resetting β")
            β = 0.0
        end

        println("β = ", β)

        force_grad_dir = false
        
        if gPrevUsed
            for i in 1:Nspin
                gPrev[i][:] = g[i][:]
            end
        end
        gKnormPrev = gKnorm

        # Update search direction
        for i in 1:Nspin
            @views d[i][:] = -Kg[i][:] + β*d[i][:]
        end

        #constrain_search_dir!( d, psi, dVol )

        α = linmin_grad!( Ham, psis, g, d, Rhoe, psisc, gt )

        # Update psi
        for i in 1:Nspin
            @views psis[i][:] = psis[i][:] + α*d[i][:]
            ortho_sqrt!(psis[i])
            psis[i] = psis[i]/sqrt(dVol)
        end

        calc_rhoe!(Ham, psis, Rhoe)
        update!(Ham, Rhoe)

        # Calc new Etot and gradient
        Etot = calc_energies_grad!( Ham, psis, Rhoe, g, Kg, Hsub )
        
        diffE = Etot_old - Etot
        @printf("Emin_PCG step %8d = %18.10f  %12.7e\n", iter, Etot, diffE)
        if diffE < 0.0
            println("*** WARNING: Etot is not decreasing")
        end

        if abs(diffE) < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        if Nconverges >= 2
            @printf("\nEmin_PCG_dot is converged in iter: %d\n", iter)
            break
        end

        Etot_old = Etot

    end

    evals = zeros(Float64,Nstates,Nspin)
    # Calculate eigenvalues
    for i in 1:Nspin
        evals[:,i], evecs = eigen(Hermitian(Hsub[i]))
        psis[i] = psis[i]*evecs
    end

    @printf("\n")
    @printf("----------------------------\n")
    @printf("Final Kohn-Sham eigenvalues:\n")
    @printf("----------------------------\n")
    @printf("\n")
    Focc = Ham.electrons.Focc
    for i in 1:Nspin
        @printf("ispin = %d\n", i)
        for ist in 1:Nstates
            @printf("%3d %18.10f %18.10f\n", ist, Focc[ist,i], evals[ist,i])
        end
    end

    @printf("\n")
    @printf("----------------------------\n")
    @printf("Final Kohn-Sham energies:\n")
    @printf("----------------------------\n")
    @printf("\n")
    println(Ham.energies, banner=false)

    return

end
