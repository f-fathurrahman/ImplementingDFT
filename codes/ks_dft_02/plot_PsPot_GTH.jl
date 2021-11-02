push!(LOAD_PATH, pwd())

using Printf

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

using KSDFT02Module

function do_plot(psp::PsPot_GTH;
    NptsPlot=500, rmax=10.0, rmax_plot=5.0, filesave="IMG_psp"
)
    
    r = range(0.0, stop=rmax, length=NptsPlot)
    
    Vloc = zeros(NptsPlot)
    for i in 1:NptsPlot
        Vloc[i] = eval_Vloc_R(psp, r[i])
    end

    plt.clf()
    plt.plot(r, Vloc)
    plt.xlim(0.0, rmax_plot)
    plt.grid()
    plt.savefig(filesave*"_Vloc.pdf")


end

#do_plot(PsPot_GTH("../pseudopotentials/pade_gth/H-q1.gth"), filesave="IMG_pade_H_q1", rmax_plot=1.0)
#do_plot(PsPot_GTH("../pseudopotentials/pade_gth/He-q2.gth"), filesave="IMG_pade_He_q2", rmax_plot=1.0)
do_plot(PsPot_GTH("../pseudopotentials/pade_gth/Li-q3.gth"), filesave="IMG_pade_Li_q3", rmax_plot=1.0)