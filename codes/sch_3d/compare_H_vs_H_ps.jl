using SpecialFunctions

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

function pot_H_coul( r::Float64 )
    return -1.0/r
end

function pot_Hps_GTH( r::Float64 )
    # Parameters
    Zval = 1
    rloc = 0.2
    C1 = -4.0663326
    C2 = 0.6678322
    rrloc = r/rloc
    V = -Zval/r * erf( r/(sqrt(2.0)*rloc) ) + (C1 + C2*rrloc^2)*exp(-0.5*(rrloc)^2)
    return V
end

function main()

    Npoints = 1000
    r = range(1e-2, stop=1.5, length=Npoints)
    V_coul = pot_H_coul.(r)
    V_H_ps = pot_Hps_GTH.(r)
    
    #println(V_coul)
    #println(V_H_ps)
    plt.clf()
    plt.plot(r, V_coul, label="Coulomb")
    plt.plot(r, V_H_ps, label="pspot")
    plt.ylim(-10.0, 0.0)
    plt.legend()
    plt.grid()
    plt.savefig("IMG_H_Coulomb_vs_pspot.pdf")

end

main()