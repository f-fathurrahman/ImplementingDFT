import PyPlot
const plt = PyPlot

function main()
    ω = 0.22
    x = range(-25.0, stop=25.0, length=51)
    V = 0.5 * ω^2 * x.^2

    plt.plot(x, V)
    plt.savefig("IMG_harm_pot.pdf")
end

main()