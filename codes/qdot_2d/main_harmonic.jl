using Printf
using LinearAlgebra
using SparseArrays
using Random
using Gnuplot

include("INC_qdot_2d.jl")

function pot_harmonic( grid::FD2dGrid; ω=1.0 )
    Npoints = grid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = grid.r[1,i]
        y = grid.r[2,i]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 )
    end
    return Vpot
end

function main()

    Random.seed!(1234)

    Nx = 81
    Ny = 81
    L = 50.0
    grid = FD2dGrid( (-L/2,L/2), Nx, (-L/2,L/2), Ny )

    ∇2 = build_nabla2_matrix( grid )

    Vpot = pot_harmonic( grid, ω=0.22 )
    
    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    prec = ILU0Preconditioner(Ham)

    @printf("sizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)
    @printf("sizeof prec = %18.10f MiB\n", Base.summarysize(prec)/1024/1024)

    dVol = grid.dVol
    Nstates = 5
    Npoints = Nx*Ny
    X = rand(Float64, Npoints, Nstates)
    ortho_sqrt!(X, dVol)

    println("Check normalization")
    for i in 1:Nstates
        @printf("dot: %18.10f\n", dot(X[:,i],X[:,i])*dVol)
    end

    println("Check orthogonal")
    ist = 1
    for i in 1:Nstates
        if i != ist
            @printf("dot: %18.10f\n", dot(X[:,i],X[:,ist])*dVol)
        end
    end

    evals = diag_Emin_PCG!( Ham, X, prec, verbose=true )
    #evals = diag_LOBPCG!( Ham, X, prec, verbose=true, tol=1e-10 )

    X = X/sqrt(grid.dVol) # renormalize

    println("Check normalization")
    for i in 1:Nstates
        @printf("dot: %18.10f\n", dot(X[:,i],X[:,i])*dVol)
    end

    println("Check orthogonal")
    ist = 1
    for i in 1:Nstates
        if i != ist
            @printf("dot: %18.10f\n", dot(X[:,i],X[:,ist])*dVol)
        end
    end

    @printf("\n\nEigenvalues\n")
    for i in 1:Nstates
        @printf("%5d %18.10f\n", i, evals[i])
    end

    for ist in 1:Nstates
        filesave = "IMG_X_"*string(ist)*".pdf"
        @gp "set term pdfcairo size 12cm,13cm fontscale 0.5" :-
        @gp :- "set output '$filesave'" :-
        @gp :- "set view 70, 40" :-
        @gsp :- grid.x grid.y reshape(X[:,ist],Nx,Ny) "w pm3d notitle"
        @gsp :- "set xrange [-9:9]" :-
        @gsp :- "set yrange [-9:9]"
    end


end

main()



