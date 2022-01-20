function op_Vtot(H::Hamiltonian, x::Array{Float64,1})
    # application of the potential part of the Hamiltonian to a vector
    return (H.Vtot).*x
end

function op_Vtot(H::Hamiltonian, X::Array{Float64,2})
    # application of the potential part of the Hamiltonian to a matric
    return broadcast(*,H.Vtot, X)
end
