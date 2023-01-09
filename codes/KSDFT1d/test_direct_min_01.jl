push!(LOAD_PATH, "./")

using Printf
using LinearAlgebra
using KSDFT1d
import Random

function create_atoms()
    Natoms = 2
    σ = ones(Float64, Natoms)*(1.0)
    masses = ones(Float64, Natoms)*42000.0
    Zvals = ones(Float64, Natoms)*2
    L = 10.0
    atpos = zeros(Float64, Natoms)
    atpos[1] = -1.0
    atpos[2] =  1.0
    return Atoms1d( atpos, Zvals, σ, masses, L )
end

function pot_gaussian( x, x0 )
    return -exp(-(x - x0)^2)/sqrt(1.0/π)^3
end

function init_Vions!(Ham)
    atpos = Ham.atoms.positions
    Natoms = Ham.atoms.Natoms
    for ia in 1:Natoms
        Ham.potentials.Ions[:] += pot_gaussian.(Ham.grid.x, atpos[ia])
    end
    return
end

function init_Hamiltonian()
    atoms = create_atoms()
    Ham = Hamiltonian1d(atoms, 51)
    init_Vions!(Ham)
    return Ham
end

# Evaluate gradient at psi
# Ham matrix is not updated
function calc_grad( Ham, psi; update=false )

    if update
        update_from_wavefunc!(Ham, psi)
    end
    
    ispin = 1
    Focc = Ham.electrons.Focc
    hx = Ham.grid.hx

    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    #
    grad = zeros(Float64, Nbasis, Nstates)

    Hpsi = Ham*psi
    for i in 1:Nstates
        @views grad[:,i] .= Hpsi[:,i]
        for j in 1:Nstates
            @views λ = dot(psi[:,j], Hpsi[:,i]) * hx
            @views grad[:,i] .= grad[:,i] .-  λ * psi[:,j]
        end
        grad[:,i] .= Focc[i,ispin]*grad[:,i]
    end

    # the usual case of constant occupation numbers
    if all(Focc[:,ispin] .== 2.0) == true || all(Focc[:,ispin] .== 1.0) == true
        # immediate return
        return grad
    end

    #F = Matrix(Diagonal(Focc[:,ispin]))
    #ℍ = psi' * Hpsi * hx
    #HFH = ℍ*F - F*ℍ
    #ℚ = 0.5*HFH
    #@views grad[:,:] .+= psi*ℚ # additional contributions
    
    return grad
end

function ortho_sqrt( psi::Array{Float64,2} )
    Udagger = inv(sqrt(psi'*psi))
    return psi*Udagger
end

function ortho_sqrt!( psi::Array{Float64,2} )
    Udagger = inv(sqrt(psi'*psi))
    psi[:,:] = psi*Udagger
    return
end

function update_from_wavefunc!(Ham, psi)
    Npoints = Ham.grid.Npoints
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    Vtot = Ham.potentials.Total
    rhoe = Ham.rhoe

    # Electron density
    calc_rhoe!(Ham, psi, rhoe)
    #println("integ rhoe = ", sum(rhoe)*hx)

    # Update the potentials
    # Note that Vxc, Vhartree, and Vtot refers to Ham.potentials
    ρ = reshape(rhoe, Npoints) # FIXME: need to do this is Poisson_solve_sum!
    Poisson_solve_sum!(Ham.grid, ρ, Vhartree)
    Vxc[:] = calc_Vxc_1d(Ham.xc_calc, rhoe)
    @views Vtot[:,1] = Vion[:] + Vhartree[:] + Vxc[:,1]

    return
end

# Ham will be modified (the potential terms)
function calc_KohnSham_Etotal!(Ham, psi)
    
    update_from_wavefunc!(Ham, psi)

    hx = Ham.grid.hx    
    Npoints = Ham.grid.Npoints
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    Vtot = Ham.potentials.Total
    rhoe = Ham.rhoe

    # Evaluate total energy components
    Ekin = calc_E_kin(Ham, psi)
    Ham.energies.Kinetic = Ekin
    
    Ehartree = 0.5*dot(rhoe[:,1], Vhartree)*hx
    Ham.energies.Hartree = Ehartree

    Eion = dot(rhoe, Vion)*hx
    Ham.energies.Ion = Eion

    epsxc = calc_epsxc_1d(Ham.xc_calc, rhoe[:,1])
    Exc = dot(rhoe, epsxc)*hx
    Ham.energies.XC = Exc

    # The total energy (also include nuclei-nuclei or ion-ion energy)
    Etot = Ekin + Ehartree + Eion + Exc + Ham.energies.NN

    return Etot
end

function prec_invK(Ham::Hamiltonian1d, v)
    return inv(Ham.Kmat)*v
end


# Ideal preconditioner (expensive to calculate)
function prec_invHam(Ham::Hamiltonian1d, v)
    Kmat = Ham.Kmat
    Vtot = Ham.potentials.Total
    Hmat = Kmat + diagm( 0 => Vtot[:,1] )
    λ = eigvals(Hmat)
    vout = similar(v)
    for i in 1:size(v,2)
        @views vout[:,i] = inv(Hmat - λ[i]*I)*v[:,i]
    end
    return vout
end


function main()

    Random.seed!(1234)

    Ham = init_Hamiltonian()

    hx = Ham.grid.hx
    Npoints = Ham.grid.Npoints
    Nelectrons = Ham.electrons.Nelectrons
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc

    Kmat = Ham.Kmat
    Vtot = Ham.potentials.Total
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    rhoe = Ham.rhoe

    Nspin = 1
    Hmat = zeros(Float64, Npoints, Npoints)
    epsxc = zeros(Float64, Npoints)
    rhoe_new = zeros(Float64, Npoints, Nspin)

    Etot = Inf
    Etot_old = Etot
    E_NN = calc_E_NN(Ham.atoms)
    Ham.energies.NN = E_NN

    psi = ortho_sqrt( rand(Float64, Npoints, Nstates) )
    psi[:] = psi[:]/sqrt(hx)
    display(psi' * psi * hx)


    α_t = 3e-5
    g = zeros(Float64, Npoints, Nstates)
    Kg = zeros(Float64, Npoints, Nstates)
    gt = zeros(Float64, Npoints, Nstates)
    d = zeros(Float64, Npoints, Nstates)
    psic = zeros(Float64, Npoints, Nstates)

    g_old = zeros(Float64, Npoints, Nstates)
    Kg_old = zeros(Float64, Npoints, Nstates)
    d_old = zeros(Float64, Npoints, Nstates)

    β = 0.0

    for iterEmin in 1:50
        
        Etot_old = Etot
        Etot = calc_KohnSham_Etotal!(Ham, psi)

        g[:,:] .= calc_grad(Ham, psi)

        dg = dot(g,g)*hx
        dEtot = abs(Etot - Etot_old)
        @printf("iterEmin = %3d Etot = %18.10f dEtot = %10.5e dg = %10.5e\n", 
            iterEmin, Etot, dEtot, dg)
        if dEtot < 1e-6
            println("Converged")
            break
        end

        Kg[:,:] .= prec_invK(Ham, g)
        #Kg[:,:] .= prec_invHam(Ham, g)

        if iterEmin >= 2
            β = real( dot(g - g_old, Kg) )/real( dot(g_old,Kg_old) )
            if β < 0.0
                β = 0.0
            end
            #println("β = ", β)
        end

        @views d[:] .= -Kg[:] + β*d_old[:]
        psic[:] = ortho_sqrt( psi + α_t*d )
        psic[:] = psic[:]/sqrt(hx)

        #_ = calc_KohnSham_Etotal!(Ham, psic)
        gt[:,:] = calc_grad(Ham, psic; update=true)

        denum = real(sum(conj(g-gt).*d))
        if denum != 0.0
            α = abs( α_t*real(sum(conj(g).*d))/denum )
        else
            α = 0.0
        end
        #denum = sum( (g - gt) .* d )
        #if denum != 0.0
        #    α = abs( α_t * sum(g .* d)/denum )
        #else
        #    α = 0.0
        #end
        #println("α = ", α)

        # Update wavefunction
        psi[:,:] .= ortho_sqrt(psi + α*d)
        psi[:] = psi[:]/sqrt(hx)

        # Previous
        @views g_old[:] .= g[:]
        @views Kg_old[:] .= Kg[:]
        @views d_old[:] .= d[:]
    end

    Hr = psi' * (Ham*psi) * hx
    band_energies = eigvals(Hr)
    display(band_energies); println()

    println(Ham.energies)

end

main()
#@time main()