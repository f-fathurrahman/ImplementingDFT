function linemin_armijo(Ham, psi, Haux, d_in, d_Haux_in, E1;
    α0=10.0, reduce_factor=0.5, α_safe=1e-8
)
    #
    hx = Ham.grid.hx
    #
    τ = reduce_factor # reduction factor
    #
    α = α0
    d = α*d_in
    d_Haux = α*d_Haux_in
    #
    psi_new = psi + d
    Haux_new = Haux + d_Haux
    #
    Udagger = inv(sqrt(psi_new'*psi_new)) ./ sqrt(hx)
    psi_new[:,:] = psi_new*Udagger
    Haux_new = Udagger' * Haux_new * Udagger
    Urot = transform_psi_Haux!(psi_new, Haux_new)
    #
    E_new = calc_Lfunc_Haux!(Ham, psi_new, Haux_new)
    #
    is_success = false
    for iterE in 1:20
        #
        dE = E_new - E1
        #@printf("LineMin: %3d %18.5e %18.5e\n", iterE, α, dE)
        if dE < 0.0
            println("E_new is smaller, α=$(α) is accepted, found in iterE=$(iterE)")
            is_success = true
            break
        end        
        #
        α = τ*α
        d[:,:] = α*d_in
        d_Haux[:,:] = α*d_Haux_in
        #
        psi_new[:,:] = psi + d
        Haux_new[:,:] = Haux + d_Haux
        #
        Udagger[:,:] = inv(sqrt(psi_new'*psi_new)) ./ sqrt(hx)
        psi_new[:,:] = psi_new*Udagger
        Haux_new[:,:] = Udagger' * Haux_new * Udagger
        Urot[:,:] = transform_psi_Haux!(psi_new, Haux_new)
        #
        E_new = calc_Lfunc_Haux!(Ham, psi_new, Haux_new)
    end

    if is_success
        return α, true
    else
        #println("LineMin: unsuccessful, returning last α")
        #return α
        #error("LineMin: unsuccessful")
        println("LineMin: unsuccessful, returning α_safe")
        return α_safe, true #false
    end
end