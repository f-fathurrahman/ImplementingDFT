function create_Ham_Al_atom( N::Int64; grid_type=:FD, Nstates_extra=4, manual_Focc=false)
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
    Ham = Hamiltonian( atoms, pspfiles, grid,
        Nstates_extra=Nstates_extra, Nspin=2 )
    if manual_Focc
        # HACK, not yet found a good API to set this
        Ham.electrons.Focc[1:2,1] = [1.0, 0.5]
        Ham.electrons.Focc[1:2,2] = [1.0, 0.5]
    end
    #
    return Ham
end


function create_Ham_Al2( N::Int64; grid_type=:FD, Nstates_extra=4)
    atoms = Atoms( xyz_string=
        """
        2

        Al   -1.0  0.0  0.0
        Al    1.0  0.0  0.0
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
    Ham = Hamiltonian( atoms, pspfiles, grid,
        Nstates_extra=Nstates_extra, Nspin=2 )
    return Ham
end

function create_Ham_C_atom( N::Int64; grid_type=:FD, Nstates_extra=4)
    atoms = Atoms( xyz_string=
        """
        1

        C   0.0  0.0  0.0
        """) # coordinates are in angstrom
    pspfiles = [ joinpath(DIR_PSP,"C-q4.gth") ]
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
    return Ham
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
