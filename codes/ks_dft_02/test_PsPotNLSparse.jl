push!(LOAD_PATH, pwd())

using Printf
using Random
using LinearAlgebra
using SpecialFunctions

using KSDFT02Module

const DIR_PSP = "../pseudopotentials/pade_gth/"

function main()
    atoms = Atoms( xyz_string=
        """
        2

        Al  1.06196 0.0  0.0
        Al -1.06196 0.0  0.0
        """)
    pspfiles = [ joinpath(DIR_PSP,"Al-q3.gth") ]
    Nspecies = atoms.Nspecies
    pspots = Array{PsPot_GTH}(undef,Nspecies)
    for isp = 1:Nspecies
        pspots[isp] = PsPot_GTH( pspfiles[isp] )
    end

    AA = -8.0*ones(3)
    BB =  8.0*ones(3)
    NN = [40,40,40]
    grid = FD3dGrid( NN, AA, BB )
    println(grid)

    pspotNL = PsPotNL( atoms, pspots, grid )
    check_betaNL_norm( grid, pspotNL )
end

main()