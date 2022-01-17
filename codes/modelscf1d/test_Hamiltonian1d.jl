push!(LOAD_PATH, "./")

using Printf
using ModelSCF1d

function init_atoms()
    Natoms = 8
    σ = ones(Float64, Natoms)*(1.0)
    masses = ones(Float64, Natoms)*42000.0
    Zvals = ones(Float64, Natoms)*2
    L = 80.0
    dx = 0.5
    atpos = zeros(Float64, Natoms)
    println("Atomic positions:")
    for ia in 1:Natoms
        atpos[ia] = (ia - 0.5)*L/Natoms + dx
        @printf("%3d %18.10f\n", ia, atpos[ia])
    end
    return Atoms1d( atpos, Zvals, σ, masses, L )
end

function main()
    atoms = init_atoms()
    
    dx_in = 0.5
    Nstates_extra = 10
    κ = 0.0100
    ε0 = 10.0
    Ham = Hamiltonian1d(atoms, dx_in, κ, ε0, Nstates_extra=Nstates_extra)
    #println(Ham)

    Hmat = get_matrix(Ham)
    println("size Hmat = ", size(Hmat))
    for ip in 1:10
        @printf("%3d %18.10f\n", ip, Hmat[ip])
    end
end

main()