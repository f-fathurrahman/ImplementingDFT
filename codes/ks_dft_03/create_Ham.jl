function create_Ham_Al_atom( N::Int64; grid_type=:FD )
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
    return Hamiltonian( atoms, pspfiles, grid, Nstates_extra=4 )
end

function create_Ham_C_atom( N::Int64; grid_type=:FD )
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
    return Hamiltonian( atoms, pspfiles, grid, Nstates_extra=4 )
end


function create_Ham_LiH( N::Int64; grid_type=:FD, Nstates_extra=0 )
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES,"LiH.xyz") )
    pspfiles = [ joinpath(DIR_PSP,"H-q1.gth"),
                 joinpath(DIR_PSP,"Li-q1.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    if (grid_type == :sinc) || (grid_type == :LF)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB )
    end
    return Hamiltonian( atoms, pspfiles, grid, Nstates_extra=Nstates_extra ) 
end

function create_Ham_CO( N::Int64; grid_type=:FD, Nstates_extra=0 )
    atoms = Atoms( xyz_string=
        """
        2

        C   0.575  0.0  0.0
        O  -0.575  0.0  0.0
        """) # coordinates are in angstrom
    pspfiles = [ joinpath(DIR_PSP,"C-q4.gth"),
                 joinpath(DIR_PSP,"O-q6.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    if (grid_type == :sinc) || (grid_type == :LF)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB )
    end
    return Hamiltonian( atoms, pspfiles, grid, Nstates_extra=Nstates_extra )  
end

function create_Ham_H2O( N::Int64; grid_type=:FD, Nstates_extra=0 )
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES,"H2O.xyz") )
    pspfiles = [ joinpath(DIR_PSP,"O-q6.gth"),
                 joinpath(DIR_PSP,"H-q1.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    if (grid_type == :sinc) || (grid_type == :LF)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB )
    end
    Ham = Hamiltonian( atoms, pspfiles, grid, Nstates_extra=Nstates_extra )
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
    Ham = Hamiltonian( atoms, pspfiles, grid, Nstates_extra=Nstates_extra )
end

function create_Ham_H( N::Int64; grid_type=:FD, Nstates_extra=0 )
    atoms = Atoms( xyz_string=
        """
        1

        H  0.0  0.0  0.0
        """ )
    pspfiles = [joinpath(DIR_PSP,"H-q1.gth")]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    if (grid_type == :sinc) || (grid_type == :LF)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB )
    end
    return Hamiltonian( atoms, pspfiles, grid, Nstates_extra=Nstates_extra )
end

function create_Ham_H2( N::Int64; grid_type=:FD, Nstates_extra=0 )
    atoms = Atoms( xyz_string=
        """
        2

        H   0.75  0.0  0.0
        H  -0.75  0.0  0.0
        """, in_bohr=true )
    pspfiles = [joinpath(DIR_PSP,"H-q1.gth")]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    if (grid_type == :sinc) || (grid_type == :LF)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB )
    end
    return Hamiltonian( atoms, pspfiles, grid )
end

function create_Ham_NH3( N::Int64; grid_type=:FD, Nstates_extra=0 )
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES,"NH3.xyz") )
    pspfiles = [ joinpath(DIR_PSP,"N-q5.gth"),
                 joinpath(DIR_PSP,"H-q1.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    if (grid_type == :sinc) || (grid_type == :LF)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB )
    end
    return Hamiltonian( atoms, pspfiles, grid, Nstates_extra=Nstates_extra )
end

function create_Ham_CH4( N::Int64; grid_type=:FD, Nstates_extra=0 )
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES,"CH4.xyz") )
    pspfiles = [ joinpath(DIR_PSP,"C-q4.gth"),
                 joinpath(DIR_PSP,"H-q1.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    if (grid_type == :sinc) || (grid_type == :LF)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB )
    end
    return Hamiltonian( atoms, pspfiles, grid, Nstates_extra=Nstates_extra )
end

function create_Ham_Ar( N::Int64; grid_type=:FD, Nstates_extra=0 )
    atoms = Atoms( xyz_string=
        """
        1

        Ar  0.0  0.0  0.0
        """ )
    pspfiles = [ joinpath(DIR_PSP,"Ar-q8.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    if (grid_type == :sinc) || (grid_type == :LF)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB )
    end
    return Hamiltonian( atoms, pspfiles, grid, Nstates_extra=Nstates_extra )
end

function create_Ham_Ne( N::Int64; grid_type=:FD, Nstates_extra=0 )
    atoms = Atoms( xyz_string=
        """
        1

        Ne  0.0  0.0  0.0
        """ )
    pspfiles = [joinpath(DIR_PSP,"Ne-q8.gth")]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    if (grid_type == :sinc) || (grid_type == :LF)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB )
    end
    return Hamiltonian( atoms, pspfiles, grid, Nstates_extra=Nstates_extra )
end

function create_Ham_SiH4( N::Int64; grid_type=:FD, Nstates_extra=0 )
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES,"SiH4.xyz") )
    pspfiles = [ joinpath(DIR_PSP,"Si-q4.gth"),
                 joinpath(DIR_PSP,"H-q1.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    grid = FD3dGrid( NN, AA, BB )
    if (grid_type == :sinc) || (grid_type == :LF)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB )
    end
    return Hamiltonian( atoms, pspfiles, grid, Nstates_extra=Nstates_extra )
end

function create_Ham_HCl( N::Int64; grid_type=:FD, Nstates_extra=0 )
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES,"HCl.xyz") )
    pspfiles = [ joinpath(DIR_PSP,"Cl-q7.gth"),
                 joinpath(DIR_PSP,"H-q1.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    grid = FD3dGrid( NN, AA, BB )
    if (grid_type == :sinc) || (grid_type == :LF)
        grid = LF3dGrid( NN, AA, BB, types=(:sinc,:sinc,:sinc) )
    else
        grid = FD3dGrid( NN, AA, BB )
    end
    return Hamiltonian( atoms, pspfiles, grid, Nstates_extra=Nstates_extra )
end
