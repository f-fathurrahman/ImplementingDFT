function create_Ham_Al_atom( N::Int64; grid_type=:FD, Nstates_extra=4)
    atoms = Atoms( xyz_string=
        """
        1

        Al   0.0  0.0  0.0
        """) # coordinates are in angstrom
    pspfiles = [ joinpath(DIR_PSP,"Al-q3.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    if (grid_type == :sinc) || (grid_type == :LF)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB )
    end
    return Hamiltonian( atoms, pspfiles, grid,
        Nstates_extra=Nstates_extra, Nspin=2 )
end

function create_Ham_O2( N::Int64; grid_type=:FD, Nstates_extra=0 )
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES,"O2.xyz") )
    pspfiles = [ joinpath(DIR_PSP,"O-q6.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    if (grid_type == :sinc) || (grid_type == :LF)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB )
    end
    Ham = Hamiltonian( atoms, pspfiles, grid,
        Nstates_extra=Nstates_extra, Nspin=2 )
end
