# Ham will be modified (the potential terms)
function calc_KohnSham_Etotal!(
    Ham,
    psi::Matrix{Float64},
    Haux::Matrix{Float64}
)

    Urot = update_from_Haux!(Ham, Haux)
    update_from_wavefunc!(Ham, psi*Urot)

    hx = Ham.grid.hx    
    Npoints = Ham.grid.Npoints
    Vion = Ham.potentials.Ions
    Vxc = Ham.potentials.XC
    Vhartree = Ham.potentials.Hartree
    Vtot = Ham.potentials.Total
    rhoe = Ham.rhoe
    mTS = Ham.energies.mTS

    # Evaluate total energy components
    Ekin = calc_E_kin(Ham, psi) # XXX: use psi*Urot ?
    Ham.energies.Kinetic = Ekin
    
    Ehartree = 0.5*dot(rhoe[:,1], Vhartree)*hx
    Ham.energies.Hartree = Ehartree

    Eion = dot(rhoe, Vion)*hx
    Ham.energies.Ion = Eion

    epsxc = calc_epsxc_1d(Ham.xc_calc, rhoe[:,1])
    Exc = dot(rhoe, epsxc)*hx
    Ham.energies.XC = Exc

    # The total energy (also include nuclei-nuclei or ion-ion energy)
    Etot = Ekin + Ehartree + Eion + Exc + Ham.energies.NN + mTS

    println(Ham.energies)

    return Etot
end