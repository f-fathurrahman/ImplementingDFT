using Printf
using LinearAlgebra
using Serialization

include("../common/constants.jl")
include("../common/Atoms.jl")
include("../common/XSF_utils.jl")

function main()

    psi = deserialize("wavefunc.data")

    AA = [-8.0, -8.0, -8.0]
    BB = [ 8.0,  8.0,  8.0]
    NN = [50, 50, 50]

    atoms = Atoms( xyz_string=
        """
        2

        H   0.75  0.0  0.0
        H  -0.75  0.0  0.0
        """, in_bohr=true )

    Lx = BB[1] - AA[1]; Ly = BB[2] - AA[2]; Lz = BB[3] - AA[3]
    LL = [ Lx 0.0 0.0;
           0.0 Ly 0.0;
           0.0 0.0 Lz]
    hx = Lx/(NN[1]-1)
    hy = Ly/(NN[2]-1)
    hz = Lz/(NN[3]-1)
    grid_origin = AA + 0.5*[hx,hy,hz]

    for i in 1:5
        filnam = "OUT_wavefunc_"*string(i)*".xsf"
        write_xsf(filnam, LL/ANG2BOHR, atoms.positions/ANG2BOHR, atsymbs=atoms.atsymbs, molecule=true)
        write_xsf_data3d_crystal(filnam, NN, LL/ANG2BOHR, psi[:,i], origin=grid_origin/ANG2BOHR)
    end

end

main()