push!(LOAD_PATH, pwd())

using Printf
using MyModule

function main()
    rho = [1.0, 2.0, 3.0, 4.0]
    
    E_xc, V_xc = calc_E_xc_V_xc_2d(rho)
    println("E_xc = ", E_xc)
    println("V_xc = ", V_xc)

    E_xc = calc_E_xc_2d(rho)
    println("E_xc = ", E_xc)

    V_xc = calc_V_xc_2d(rho)
    println("V_xc = ", V_xc)

end

main()