push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using MyModule

function pot_harmonic( grid; ω=1.0, center=[0.0, 0.0] )
    Npoints = grid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = grid.r[1,i] - center[1]
        y = grid.r[2,i] - center[2]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 )
    end
    return Vpot
end

include("KS_solve_Emin_PCG.jl")

function main()

    Random.seed!(1234)

    AA = [-25.0, -25.0]
    BB = [ 25.0,  25.0]
    NN = [81, 81]

    grid = FD2dGrid( NN, AA, BB )
    #grid = LF2dGrid( NN, AA, BB, types=(:sinc,:sinc) )

    V_ext = pot_harmonic( grid, ω=0.22 )
    
    Nstates = 1
    Nelectrons = 2*Nstates
    Ham = Hamiltonian( grid, V_ext, Nelectrons=Nelectrons )

    @printf("sizeof Ham  = %18.10f MiB\n", Base.summarysize(Ham)/1024/1024)

    Nbasis = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol

    psi = rand(Float64,Nbasis,Nstates)
    ortho_sqrt!(psi,dVol)

    KS_solve_Emin_PCG!(Ham, psi)

end

@time main()
