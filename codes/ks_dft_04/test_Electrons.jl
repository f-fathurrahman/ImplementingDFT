push!(LOAD_PATH, pwd())

using Printf

using MyModule

const DIR_PSP = "../pseudopotentials/pade_gth/"
const DIR_STRUCTURES = "../structures"

function test_O2()

    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES,"O2.xyz") )
    pspfiles = [ joinpath(DIR_PSP,"O-q6.gth") ]

    Nspecies = atoms.Nspecies
    pspots = Array{PsPot_GTH}(undef,Nspecies)
    for isp = 1:Nspecies
        pspots[isp] = PsPot_GTH(pspfiles[isp])
    end

    #electrons = Electrons( atoms, pspots, (5,7) )
    #println(electrons)
    #electrons = Electrons( atoms, pspots, (7,5), Nstates_extra=1 )
    #println(electrons)

    #electrons = Electrons(atoms, pspots, 2, Nstates_extra=1)
    #println(electrons)

    electrons = Electrons(atoms, pspots, 2, Nstates_extra=1, Nspin=2)
    println(electrons)

end

test_O2()