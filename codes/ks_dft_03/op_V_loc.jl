# Local potential operator, Vloc is the local potential to be applied instead
# of potential in Ham
function op_V_loc( V_loc::Vector{Float64}, psi )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    Vpsi = zeros(Float64,Nbasis,Nstates)
    for ist in 1:Nstates, ip in 1:Nbasis
        Vpsi[ip,ist] = V_loc[ip]*psi[ip,ist]
    end
    return Vpsi
end

# Local pspot operator
function op_V_Ps_loc( Ham::Hamiltonian, psi )
    Nbasis = size(psi,1)
    Nstates = size(psi,2)
    Vpsi = zeros(Float64,Nbasis,Nstates)
    for ist in 1:Nstates, ip in 1:Nbasis
        Vpsi[ip,ist] = Ham.V_Ps_loc[ip]*psi[ip,ist]
    end
    return Vpsi
end