push!(LOAD_PATH, "./")

using Printf
using LinearAlgebra
using ModelSCF1d

function init_atoms()
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

function calc_rhoe!(
    Ham::Hamiltonian1d,
    psi::Matrix{Float64},
    rhoe::Matrix{Float64}
)
    Npoints = Ham.Ns
    Focc = Ham.electrons.Focc
    Nstates = Ham.electrons.Nstates
    fill!(rhoe, 0.0)
    ispin = 1
    for ist in 1:Nstates
        if Focc[ist,ispin] > 0.0
            for ip in 1:Npoints
                rhoe[ip,ispin] += Focc[ist,ispin] * psi[ip,ist] * psi[ip,ist]
            end
        end
    end
    return
end



function main()
    atoms = init_atoms()
    
    dx_in = 0.5
    Nstates_extra = 18
    κ = 0.0100
    ε0 = 10.0
    Ham = Hamiltonian1d(atoms, dx_in, κ, ε0, Nstates_extra=Nstates_extra)
    kT = 1/3157.7466999999997  # 1/Tbeta
    println("kT = ", kT)

    Hmat = zeros(Float64, Ham.Ns, Ham.Ns)
    Nstates = Ham.electrons.Nstates
    Nelectrons = Ham.electrons.Nelectrons

    Nspin = 1
    Focc = Ham.electrons.Focc
    ebands = Ham.electrons.ebands
    psi = zeros(Float64, Ham.Ns, Nstates) # not spin polarized
    Rhoe = Ham.rhoe
    Rhoe_new = zeros(Float64, Ham.Ns, Nspin)
    β_mix = 0.1

    for iter_scf in 1:200
        # Using dense matrix
        Hmat = get_matrix(Ham)
        evals_all, psi_all = eigen(Hmat)
        ebands[:,1] = evals_all[1:Nstates]
        psi[:,:] = psi_all[:,1:Nstates]/sqrt(Ham.dx) # normalization
        #
        #println("Check ortho psi: 1 1: ", dot(psi[:,1], psi[:,1])*Ham.dx)
        #println("Check ortho psi: 1 2: ", dot(psi[:,1], psi[:,2])*Ham.dx)
        #
        E_fermi, mTS = update_Focc!(
            Focc, smear_fermi, smear_fermi_entropy,
            ebands, Nelectrons, kT )
        #
        #for ist in 1:Nstates
        #    @printf("%3d %18.10f %18.10f\n", ist, ebands[ist,1], Focc[ist,1])
        #end
        #println("E_fermi = ", E_fermi)
        #println("mTS     = ", mTS)
        #
        calc_rhoe!(Ham, psi, Rhoe_new)
        #println("integ Rhoe_new: ", sum(Rhoe_new)*Ham.dx)
        rmseRhoe = sqrt(norm(Rhoe - Rhoe_new)/Ham.Ns)
        maeRhoe = sum(abs.(Rhoe - Rhoe_new))/Ham.Ns
        #
        # relative error
        #
        rel_err = norm(Rhoe_new - Rhoe)/norm(Rhoe)
        #
        @printf("%3d %10.5e %10.5e %10.5e\n", iter_scf, rmseRhoe, maeRhoe, rel_err)
        #
        if maeRhoe < 1e-7 && rel_err < 1e-7
            @printf("Converged\n")
            break
        end

        # Mix rhoe
        @views Rhoe[:] = β_mix*Rhoe_new[:] + (1 - β_mix)*Rhoe[:]
        #println("integ Rhoe: ", sum(Rhoe)*Ham.dx)
        #println("integ Ham.rhoe: ", sum(Ham.rhoe)*Ham.dx)
        #
        # New potential
        #
        rho_in = dropdims( sum(Rhoe, dims=2), dims=2 ) + Ham.rhoa
        @views PoissonYukawa1d_solve!( Ham.κ, Ham.ε0,
            Ham.gvec, rho_in, Ham.VHartree
        )
        Ham.Vtotal[:] = Ham.VHartree[:] # no XC potential
    end

end

main()