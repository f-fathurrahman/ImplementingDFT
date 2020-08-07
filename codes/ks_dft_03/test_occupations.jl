using SpecialFunctions: erfc
using Printf

include("smearing.jl")
include("occupations.jl")

function test_main()

    evals = [-1.2, -1.1, -1.0, -0.8, -0.7]
    Nelectrons = 6.0
    kT = 0.01

    E_fermi = find_E_fermi(smear_cold, evals, Nelectrons, kT, verbose=true)
    #println("smear_fermi = ", smear_fermi(evals[1], -1.3, kT))

    println("E_fermi = ", E_fermi)
    println("typeof(smear_fermi) = ", smear_fermi)
end

test_main()