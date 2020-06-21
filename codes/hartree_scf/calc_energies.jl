function calc_E_NN( Zvals::Array{Float64,1}, r::Array{Float64,2} )
    Natoms = length(Zvals)
    @assert Natoms == size(r,2)
    E_NN = 0.0
    for ia in 1:Natoms
        for ja in ia+1:Natoms
            dx = r[1,ja] - r[1,ia]
            dy = r[2,ja] - r[2,ia]
            dz = r[3,ja] - r[3,ia]
            r_ij = sqrt(dx^2 + dy^2 + dz^2)
            E_NN = E_NN + Zvals[ia]*Zvals[ja]/r_ij
        end
    end
    return E_NN
end

function calc_E_kin( Ham, psi::Array{Float64,2} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    E_kin = 0.0
    nabla2psi = zeros(Float64,Nbasis)
    dVol = Ham.grid.dVol
    for ist in 1:Nstates
        @views nabla2psi = -0.5*Ham.Laplacian*psi[:,ist]
        @views E_kin = E_kin + Ham.electrons.Focc[ist]*dot( psi[:,ist], nabla2psi[:] )*dVol
    end
    return E_kin
end

function calc_energies!( Ham::Hamiltonian, psi::Array{Float64,2} )
    dVol = Ham.grid.dVol
    Ham.energies.Kinetic = calc_E_kin( Ham, psi )
    Ham.energies.Ps_loc = sum( Ham.V_Ps_loc .* Ham.rhoe )*dVol
    Ham.energies.Hartree = 0.5*sum( Ham.V_Hartree .* Ham.rhoe )*dVol
    return
end

function calc_energies(Ham, psi)
    calc_energies!(Ham, psi)
    return Ham.energies
end