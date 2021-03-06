function create_Ham_H_periodic( N::Int64; grid_type=:FD )
    atoms = Atoms( xyz_string=
        """
        1

        H  8.0  8.0  8.0
        """,
        in_bohr=true,
        pbc=(true,true,true),
        LatVecs=16.0*diagm(ones(3)) )
    println(atoms)
    pspfiles = [joinpath(DIR_PSP,"H-q1.gth")]
    AA = zeros(3)
    BB = 16.0*ones(3)
    NN = [N,N,N]
    if grid_type == :LF
        grid = LF3dGrid( NN, AA, BB, types=(:P,:P,:P) )
    else
        grid = FD3dGrid( NN, AA, BB, pbc=(true,true,true) )
    end
    return Hamiltonian( atoms, pspfiles, grid )
end

function create_Ham_He_periodic( N::Int64; grid_type=:FD )
    atoms = Atoms( xyz_string=
        """
        1

        He  8.0  8.0  8.0
        """,
        in_bohr=true,
        pbc=(true,true,true),
        LatVecs=16.0*diagm(ones(3)) )
    println(atoms)
    pspfiles = [joinpath(DIR_PSP,"He-q2.gth")]
    AA = zeros(3)
    BB = 16.0*ones(3)
    NN = [N,N,N]
    if grid_type == :LF
        grid = LF3dGrid( NN, AA, BB, types=(:P,:P,:P) )
    else
        grid = FD3dGrid( NN, AA, BB, pbc=(true,true,true) )
    end
    return Hamiltonian( atoms, pspfiles, grid )
end


function create_Ham_Ne_periodic( N::Int64; grid_type=:FD )
    atoms = Atoms( xyz_string=
        """
        1

        Ne  8.0  8.0  8.0
        """,
        pbc=(true,true,true), in_bohr=true,
        LatVecs=16.0*diagm(ones(3))
    )
    pspfiles = [joinpath(DIR_PSP,"Ne-q8.gth")]
    AA = zeros(3)
    BB = 16.0*ones(3)
    NN = [N,N,N]
    if grid_type == :LF
        grid = LF3dGrid( NN, AA, BB, types=(:P,:P,:P) )
    else
        grid = FD3dGrid( NN, AA, BB, pbc=(true,true,true) )
    end
    return Hamiltonian( atoms, pspfiles, grid )
end

# This will not work correctly as the coordinate of LiH is outside the box
function create_Ham_LiH_periodic_bad( N::Int64; grid_type=:FD )
    atoms = Atoms(
        xyz_file=joinpath(DIR_STRUCTURES,"LiH.xyz"),
        pbc=(true,true,true),
        LatVecs=16.0*diagm(ones(3))
    )
    pspfiles = [ joinpath(DIR_PSP,"H-q1.gth"),
                 joinpath(DIR_PSP,"Li-q1.gth") ]
    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [N,N,N]
    if grid_type == :LF
        grid = LF3dGrid( NN, AA, BB, types=(:P,:P,:P) )
    else
        grid = FD3dGrid( NN, AA, BB, pbc=(true,true,true) )
    end
    return Hamiltonian( atoms, pspfiles, grid )    
end

function create_Ham_LiH_periodic( N::Int64; grid_type=:FD )
    atoms = Atoms(
        xyz_file=joinpath(DIR_STRUCTURES,"LiH_cell.xyz"),
        pbc=(true,true,true),
        LatVecs=16.0*diagm(ones(3))
    )
    pspfiles = [ joinpath(DIR_PSP,"H-q1.gth"),
                 joinpath(DIR_PSP,"Li-q3.gth") ]
    AA = zeros(3)
    BB = 16.0*ones(3)
    NN = [N,N,N]
    if grid_type == :LF
        grid = LF3dGrid( NN, AA, BB, types=(:P,:P,:P) )
    else
        grid = FD3dGrid( 16.0*ones(3), NN, pbc=(true,true,true) )
    end
    return Hamiltonian( atoms, pspfiles, grid )    
end

# With nonlocal psp
function create_Ham_LiH_periodic_v2( N::Int64; grid_type=:FD )
    atoms = Atoms(
        xyz_file=joinpath(DIR_STRUCTURES,"LiH_cell.xyz"),
        pbc=(true,true,true),
        LatVecs=16.0*diagm(ones(3))
    )
    pspfiles = [ joinpath(DIR_PSP,"H-q1.gth"),
                 joinpath(DIR_PSP,"Li-q1.gth") ]
    AA = zeros(3)
    BB = 16.0*ones(3)
    NN = [N,N,N]
    if grid_type == :LF
        grid = LF3dGrid( NN, AA, BB, types=(:P,:P,:P) )
    else
        grid = FD3dGrid( 16.0*ones(3), NN, pbc=(true,true,true) )
    end
    return Hamiltonian( atoms, pspfiles, grid )    
end


function create_Ham_CH4_periodic( N::Int64; grid_type=:FD )
    atoms = Atoms(
        xyz_file=joinpath(DIR_STRUCTURES,"CH4_cell.xyz"),
        pbc=(true,true,true),
        LatVecs=16.0*diagm(ones(3))
    )
    pspfiles = [ joinpath(DIR_PSP,"C-q4.gth"),
                 joinpath(DIR_PSP,"H-q1.gth") ]
    AA = zeros(3)
    BB = 16.0*ones(3)
    NN = [N,N,N]
    if grid_type == :LF
        grid = LF3dGrid( NN, AA, BB, types=(:P,:P,:P) )
    else
        grid = FD3dGrid( 16.0*ones(3), NN, pbc=(true,true,true) )
    end
    return Hamiltonian( atoms, pspfiles, grid )    
end