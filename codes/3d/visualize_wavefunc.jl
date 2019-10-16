using Printf
using LinearAlgebra
using Serialization

mutable struct Atoms
end

include("../XSF_utils.jl")

function main()

    psi = deserialize("wavefunc.data")

    LL = [10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0]
    Ns = (50, 50, 50)

    for i in 1:10
        filnam = "wavefunc_"*string(i)*".xsf"
        write_xsf(filnam, LL)
        write_xsf_data3d_crystal(filnam, Ns, LL, psi[:,i])
    end

end

main()