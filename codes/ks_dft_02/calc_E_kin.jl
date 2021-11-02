function calc_E_kin( Ham, psi::Array{Float64,2} )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    E_kin = 0.0
    dVol = Ham.grid.dVol
    nabla2psi = zeros(Float64,Nbasis)
    Focc = Ham.electrons.Focc
    for ist in 1:Nstates
        @views Kpsi = -0.5*Ham.∇2*psi[:,ist]
        @views E_kin = E_kin + Focc[ist]*dot( psi[:,ist], Kpsi[:] )*dVol
    end
    return E_kin
end