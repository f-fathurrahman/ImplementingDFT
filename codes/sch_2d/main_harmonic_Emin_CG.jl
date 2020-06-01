using Printf
using LinearAlgebra
using SparseArrays
using Random

using IncompleteLU
using AlgebraicMultigrid

import PyPlot
const plt = PyPlot

include("INC_sch_2d.jl")
include("../common/ILU0Preconditioner.jl")
include("../common/NoPreconditioner.jl")

function pot_harmonic( grid; ω=1.0 )
    Npoints = grid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = grid.r[1,i]
        y = grid.r[2,i]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 )
    end
    return Vpot
end


function check_ortho(X, dVol)
    Nstates = size(X,2)
    println("Check normalization:")
    for i in 1:Nstates
        @printf("dot: %18.10f\n", dot(X[:,i],X[:,i])*dVol)
    end

    println("Check orthogonality wrt 1st state:")
    ist = 1
    for i in 1:Nstates
        if i != ist
            @printf("dot: %18.10f\n", dot(X[:,i],X[:,ist])*dVol)
        end
    end

end

function calc_grad_evals!( Ham, ψ, g, Hsub )
    Nstates = size(ψ,2)
    Hψ = Ham*ψ
    Hsub[:] = ψ' * Hψ
    g[:,:] = Hψ - ψ*Hsub
    return
end


function main()

    Random.seed!(1234)

    Nx = 50
    Ny = 50
    grid = FD2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )

    ∇2 = build_nabla2_matrix( grid, stencil_order=9 )

    Vpot = pot_harmonic( grid )
    
    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    #prec = ilu(-0.5*∇2)
    #prec = ilu(Ham) # this should result in faster convergence
    prec = ILU0Preconditioner(Ham)
    #prec = NoPreconditioner()

    @printf("sizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)
    @printf("sizeof prec = %18.10f MiB\n", Base.summarysize(prec)/1024/1024)

    dVol = grid.dVol
    Nstates = 10
    Npoints = Nx*Ny
    X = rand(Float64, Npoints, Nstates)
    ortho_sqrt!(X, dVol)
    #check_ortho(X, dVol)

    Nbasis  = size(X,1)
    Nstates = size(X,2)

    Nstates_conv = Nstates
    @assert Nstates_conv <= Nstates_conv
    NiterMax = 1000
    α_t = 3e-5
    tol = 1e-6

    #
    # Variabls for PCG
    #
    Hsub = zeros(Float64,Nstates,Nstates)
    g = zeros(Float64,Nbasis,Nstates)
    d = zeros(Float64,Nbasis,Nstates)
    g_old = copy(d)
    d_old = copy(d)
    Kg = copy(d)
    Kg_old = copy(d)
    Xc = copy(d)
    gt = copy(d)
    
    β = 0.0
    
    Ebands_old = 0.0

    Hr = Hermitian( X' * (Ham*X) )
    evals = eigvals(Hr)
    Ebands = sum(evals)
    
    evals_old = copy(evals)
    devals = ones(Nstates)

    IS_CONVERGED = false

    for iter = 1:NiterMax

        calc_grad_evals!( Ham, X, g, Hsub )
        
        Kg[:] = g[:] # copy
        for i in 1:Nstates
            @views ldiv!(prec, Kg[:,i])
        end

        if iter != 1
            β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g_old).*Kg_old))
        end

        if β < 0.0
            β = 0.0
        end

        d = -Kg + β * d_old

        Xc = ortho_sqrt( X + α_t*d )
        calc_grad_evals!( Ham, Xc, gt, Hsub )

        denum = real(sum(conj(g-gt).*d))
        if denum != 0.0
            α = abs( α_t*real(sum(conj(g).*d))/denum )
        else
            α = 0.0
        end

        # Update wavefunction
        X[:,:] = X + α*d
        ortho_sqrt!(X)

        Hr = Hermitian( X' * ( Ham*X ) )
        evals = eigvals(Hr)
        Ebands = sum(evals)

        devals = abs.( evals - evals_old )
        evals_old = copy(evals)
        
        nconv = length( findall( devals .< tol ) )

        diffE = abs(Ebands-Ebands_old)

        @printf("CG step %8d = %18.10f   %10.7e  nconv = %5d\n", iter, Ebands, diffE, nconv)
        if nconv >= Nstates_conv
            IS_CONVERGED = true
            @printf("Emin_CG convergence: nconv = %5d in %5d iterations\n", nconv, iter)
            break
        end

        g_old = copy(g)
        d_old = copy(d)
        Kg_old = copy(Kg)
        Ebands_old = Ebands

        flush(stdout)

    end

    if !IS_CONVERGED
        @printf("\nWARNING: Emin_CG is not converged after %d iterations\n", NiterMax)
    end

    ortho_sqrt!(X)
    Hr = Hermitian( X' * (Ham*X) )
    evals, evecs = eigen(Hr)
    X[:,:] = X*evecs

    @printf("\nEigenvalues:\n\n")
    for ist = 1:Nstates
        @printf("evals[%3d] = %18.10f devals = %18.10e\n", ist, evals[ist], devals[ist] )
    end


    X = X/sqrt(grid.dVol) # renormalize

end

main()



