using SpecialFunctions: erfc, erf

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

function smear_cold_v2(x)
    MAXARG = 200.0
    xp = x - 1.0/sqrt(2.0)
    arg = min(MAXARG, xp^2)
    return 0.5 * erf(xp) + 1.0/sqrt(2.0*pi) * exp(-arg) + 0.5
end

function do_plot()

    # x = (ϵ - E_f)/(2.0*kT)
    
    smear_fermi(x) = 0.5*(1.0 - tanh(x))
    smear_gauss(x) = 0.5*erfc(x)
    smear_MP1(x) = 0.5*erfc(x) - x*exp(-x*x)/(2*sqrt(pi))
    smear_cold(x) = 0.5*erfc(x+sqrt(0.5)) + exp( -(x + sqrt(0.5))^2 )/sqrt(2*pi)

    x = range(-10.0, stop=10.0, length=500)
    plt.clf()
    plt.plot(x, smear_fermi.(x), label="Fermi")
    plt.plot(x, smear_gauss.(x), label="Gauss")
    plt.plot(x, smear_MP1.(x), label="MP1")
    plt.plot(x, smear_cold.(x), label="cold")
    plt.plot(x, smear_cold_v2.(x), label="cold v2")
    plt.grid()
    plt.legend()
    plt.savefig("IMG_smearing_func.pdf")
end
do_plot()


function smear_fermi_TS(x)
    f = 0.5*( 1 - tanh(x) )
    S = 0.0
    if f > 1e-300
        S = S - f*log(f)
    end
    if (1-f) > 1e-300
        S = S - (1-f)*log(1-f)
    end
    return S
end

function smear_cold_v2_TS(x)
    xp = x - 1.0/sqrt(2.0)
    arg = min(200.0, xp^2)
    return 1.0/sqrt(2.0*pi)*xp*exp(-arg)
end

function do_plot_entropy()

    # x = (ϵ - E_f)/(2.0*kT)
    
    smear_gauss_TS(x) = exp(-x^2) / sqrt(pi)
    smear_MP1_TS(x) = (0.5 - x*x)*exp(-x*x)/sqrt(pi)
    smear_cold_TS(x) = exp( -(x+sqrt(0.5))^2 ) * (1 + x*sqrt(2))/sqrt(pi)

    x = range(-10.0, stop=10.0, length=500)
    plt.clf()
    plt.plot(x, smear_fermi_TS.(x), label="Fermi")
    plt.plot(x, smear_gauss_TS.(x), label="Gauss")
    plt.plot(x, smear_MP1_TS.(x), label="MP1")
    plt.plot(x, smear_cold_TS.(x), label="cold")
    plt.plot(x, smear_cold_v2_TS.(x), label="cold v2")
    plt.grid()
    plt.legend()
    plt.savefig("IMG_smearing_func_entropy.pdf")
end
do_plot_entropy()