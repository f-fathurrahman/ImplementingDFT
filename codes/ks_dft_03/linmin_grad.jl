function linmin_grad_v1!(
    Ham::Hamiltonian,
    evars_::ElecVars, g::ElecGradient, d::ElecGradient, kT::Float64,
    subrot_; αt = 3e-5
)

    dVol = Ham.grid.dVol

    # Probably preallocate these outside linmin_grad?
    evars = deepcopy(evars_)
    subrot = deepcopy(subrot_)
    gt = deepcopy(g)
    Kgt = deepcopy(g)

    do_step!( Ham, αt, αt, evars, d, subrot )
    Etot = compute!( Ham, evars, gt, Kgt, kT, subrot )

    denum = dot_ElecGradient(g - gt, d, dVol)
    if denum != 0.0
        α = abs( αt * dot_ElecGradient(g, d, dVol) / denum )
    else
        α = 0.0
    end
    return α

end


function linmin_grad_v2!(
    Ham::Hamiltonian,
    evars_::ElecVars, g::ElecGradient, d::ElecGradient, kT::Float64,
    subrot_; αt = 3e-5
)

    dVol = Ham.grid.dVol

    # Probably preallocate these outside linmin_grad?
    evars = deepcopy(evars_)
    subrot = deepcopy(subrot_)
    gt = deepcopy(g)
    Kgt = deepcopy(g)

    do_step!( Ham, αt, αt, evars, d, subrot )
    Etot = compute!( Ham, evars, gt, Kgt, kT, subrot )

    denum, denum_aux = dot_ElecGradient_v2(g - gt, d, dVol)
    num, num_aux = dot_ElecGradient_v2(g, d, dVol)
    if denum != 0.0
        α = abs( αt * num/denum )
    else
        α = 0.0
    end

    if denum_aux != 0.0
        α_aux = abs( αt * num_aux/denum_aux )
    else
        α_aux = 0.0
    end

    return α, α_aux

end
