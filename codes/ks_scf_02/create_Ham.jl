function create_Ham_LiH( N::Int64 )
    atoms = Atoms( xyz_file="LiH.xyz" )
    pspfiles = [ joinpath("H-q1.gth"),
                 joinpath("Li-q1.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [41, 41, 41]
    #grid = FD3dGrid( NN, AA, BB )
    #grid = LF3dGrid( NN, AA, BB )
    grid = LF3dGrid( NN, AA, BB, type_x=:sinc, type_y=:sinc, type_z=:sinc)
    return Hamiltonian( atoms, pspfiles, grid )    
end

function create_Ham_H2O( N::Int64 )
    atoms = Atoms( xyz_file="H2O.xyz" )
    pspfiles = [ joinpath(DIR_PSP,"O-q6.gth"),
                 joinpath(DIR_PSP, "H-q1.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    grid = FD3dGrid( NN, AA, BB )
    Ham = Hamiltonian( atoms, pspfiles, grid )
end

function create_Ham_H( N::Int64 )
    atoms = Atoms( xyz_string=
        """
        1

        H  0.0  0.0  0.0
        """ )
    pspfiles = ["H-q1.gth"]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    grid = FD3dGrid( NN, AA, BB )
    return Hamiltonian( atoms, pspfiles, grid )
end

function create_Ham_NH3( N::Int64 )
    atoms = Atoms( xyz_file="NH3.xyz" )
    pspfiles = ["N-q5.gth", "H-q1.gth"]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    grid = FD3dGrid( NN, AA, BB )

    return Hamiltonian( atoms, pspfiles, grid )
end

function create_Ham_CH4( N::Int64 )
    atoms = Atoms( xyz_file="CH4.xyz" )
    pspfiles = ["C-q4.gth", "H-q1.gth"]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    grid = FD3dGrid( NN, AA, BB )
    return Hamiltonian( atoms, pspfiles, grid )
end

function create_Ham_Ar( N::Int64 )
    atoms = Atoms( xyz_string=
        """
        1

        Ar  0.0  0.0  0.0
        """ )
    pspfiles = [ joinpath(DIR_PSP, "Ar-q8.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    grid = FD3dGrid( NN, AA, BB )

    return Hamiltonian( atoms, pspfiles, grid )
end

function create_Ham_Ne( N::Int64 )
    atoms = Atoms( xyz_string=
        """
        1

        Ne  0.0  0.0  0.0
        """ )
    pspfiles = ["Ne-q8.gth"]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    grid = FD3dGrid( NN, AA, BB )

    return Hamiltonian( atoms, pspfiles, grid )
end

function create_Ham_SiH4(N::Int64)
    atoms = Atoms( xyz_file="SiH4.xyz" )
    pspfiles = ["Si-q4.gth", "H-q1.gth"]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    grid = FD3dGrid( NN, AA, BB )

    return Hamiltonian( atoms, pspfiles, grid )
end