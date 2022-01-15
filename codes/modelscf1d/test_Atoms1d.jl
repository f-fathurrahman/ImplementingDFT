push!(LOAD_PATH, "./")

using Printf
using ModelSCF1d

function main()
    Natoms = 8
    σ = ones(Float64, Natoms)*(1.0)
    masses = ones(Float64, Natoms)*42000.0
    Zvals = ones(Float64, Natoms)*2
    L = 80.0
    dx = 0.5
    atpos = zeros(Float64, Natoms)
    for ia in 1:Natoms
        atpos[ia] = (ia - 0.5)*L/Natoms + dx
        @printf("%3d %18.10f\n", ia, atpos[ia])
    end
    
    atoms = Atoms1d( atpos, Zvals, σ, masses, L )
end

main()