using SpecialFunctions: erfc
using Printf

include("smearing.jl")
include("occupations.jl")

function test_v1()

    evals = [-1.35, -1.15, -1.0, -0.8, -0.7]
    Nstates = size(evals,1)
    Nelectrons = 4.0
    kT = 0.02

    wks = [2.0]

    E_f = find_E_fermi(wgauss, evals, Nelectrons, kT)
    println("E_fermi = ", E_f)
    mTS = 0.0
    println("\nOriginal")
    for ist = 1:Nstates
        xx = wks[1]*kT*w1gauss( evals[ist], E_f, kT )
        mTS = mTS + xx
        @printf("xx = %18.10f\n", xx)
    end
    println("mTS = ", mTS)

    E_f = find_E_fermi(smear_fermi, evals, Nelectrons, kT)
    println("E_fermi = ", E_f)
    mTS = 0.0
    println("\nNew")
    for ist = 1:Nstates
        xx = -wks[1]*kT*smear_fermi_entropy( evals[ist], E_f, kT )
        mTS = mTS + xx
        @printf("xx = %18.10f\n", xx)
    end
    println("mTS = ", mTS)

end
test_v1()

function test_v2()

    evals = [-1.35, -1.15, -1.0, -0.8, -0.7]
    Nstates = size(evals,1)
    Nelectrons = 4.0
    kT = 0.02

    wks = [2.0]

    Focc = zeros(Nstates)
    E_f, mTS = update_Focc!( Focc, smear_fermi, smear_fermi_entropy,
                             evals, Nelectrons, kT )
    println("E_f  = ", E_f)
    println("mTS  = ", mTS)
    println("Focc = ", Focc)
end
test_v2()