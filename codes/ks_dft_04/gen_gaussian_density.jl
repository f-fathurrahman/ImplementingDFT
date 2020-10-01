function gen_gaussian_density!(
    grid, atoms::Atoms, pspots::Vector{PsPot_GTH}, Rhoe::Array{Float64,1}
)
    Natoms = atoms.Natoms
    Npoints = grid.Npoints

    Zvals = get_Zvals(pspots)

    σ = 0.5
    nrmfct = (2*pi*σ^2)^1.5
    Rhoe .= 0.0
    for ia in 1:Natoms
        isp = atoms.atm2species[ia]
        for ip in 1:Npoints
            dr1 = grid.r[1,ip] - atoms.positions[1,ia]
            dr2 = grid.r[2,ip] - atoms.positions[2,ia]
            dr3 = grid.r[3,ip] - atoms.positions[3,ia]
            r = sqrt(dr1^2 + dr2^2 + dr3^2)
            Rhoe[ip] = Rhoe[ip] + Zvals[isp] * exp( -r^2 / (2.0*σ^2) ) / nrmfct
        end
    end
    println("Gaussian Rhoe integral = ", sum(Rhoe)*grid.dVol)
    return
end