push!(LOAD_PATH, pwd())

using Printf

using KSDFT02Module

function main()
    psp1 = PsPot_GTH("../pseudopotentials/pade_gth/Ar-q8.gth")
    psp2 = PsPot_GTH_octopus("../pseudopotentials/pade_gth/Ar-q8.gth")

    println(psp1)
    println(psp2)
end

main()