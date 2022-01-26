push!(LOAD_PATH, "./")

using Printf
using KSDFT1d

function create_atoms()
    Natoms = 8
    σ = ones(Float64, Natoms)*(1.0)
    masses = ones(Float64, Natoms)*42000.0
    Zvals = ones(Float64, Natoms)*2
    L = 80.0
    dx = 0.5
    atpos = zeros(Float64, Natoms)
    for ia in 1:Natoms
        atpos[ia] = (ia - 0.5)*L/Natoms + dx
    end
    return Atoms1d( atpos, Zvals, σ, masses, L )
end

function main()
    atoms = create_atoms()
    Ham = Hamiltonian1d(atoms, 51)
end

main()
