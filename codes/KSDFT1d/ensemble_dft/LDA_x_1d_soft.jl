using QuadGK: quadgk
using SpecialFunctions: besselk


function FT_inter(x)
    return 2.0*besselk(0, x);
end

func1(x) = FT_inter(x) # need this?
func2(x) = x*FT_inter(x)


function my_piecewise5(c1, x1, c2, x2, x3)
    #( (c1) ? (x1) : ( (c2) ? (x2) : (x3)) )
    if c1
        return x1
    elseif c2
        return x2
    else
        return x3
    end
end


function my_piecewise3(c, x1, x2)
    #((c) ? (x1) : (x2))
    if c
        return x1
    else
        return x2
    end
end


function calc_eps_x_1d_soft(ρ)

    params_beta = 1.0
    # This is probably enough
    p_dens_threshold = 10*eps()
    p_zeta_threshold = 10*eps()

    t3 = 0.1e1 <= p_zeta_threshold;
    t4 = ρ / 0.2e1 <= p_dens_threshold || t3;
    t5 = p_zeta_threshold - 0.1e1;
    
    t7 = my_piecewise5(t3, t5, t3, -t5, 0);
    t8 = 0.1e1 + t7;

    t11 = t8 * π * params_beta * ρ;
    
    if t11 > eps()
        t12, _ = quadgk(func1, 0.0, t11, rtol=eps());
        t14, _ = quadgk(func2, 0.0, t11, rtol=eps());
    else
        # Is this OK?
        return 0.0
    end

    t15 = 0.1e1 / π;
    t16 = t14 * t15;
    t17 = 0.1e1 / params_beta;
    t18 = 0.1e1 / ρ;
    t19 = t17 * t18;
    t24 = my_piecewise3(t4, 0, -0.79577471545947667883e-1 * (t8 * t12 - t16 * t19) * t17);
    tzk0 = 0.2e1 * t24;

    return tzk0
end


function calc_V_x_1d_soft( ρ::Float64 )

    params_beta = 1.0
    p_dens_threshold = 10*eps()
    p_zeta_threshold = 10*eps()

    t3 = 0.1e1 <= p_zeta_threshold;
    t4 = ρ / 0.2e1 <= p_dens_threshold || t3;
    t5 = p_zeta_threshold - 0.1e1;
    t7 = my_piecewise5(t3, t5, t3, -t5, 0);
    t8 = 0.1e1 + t7;
    
    t11 = t8 * π * params_beta * ρ;
    if t11 > eps()
        t12, _ = quadgk(func1, 0.0, t11, rtol=eps());
        t14, _ = quadgk(func2, 0.0, t11, rtol=eps());
    else
        # Is this OK?
        return 0.0
    end

    t15 = 0.1e1 / π;
    t16 = t14 * t15;
    t17 = 0.1e1 / params_beta;
    t18 = 0.1e1 / ρ;
    t19 = t17 * t18;
    
    t24 = my_piecewise3(t4, 0, -0.79577471545947667883e-1 * (t8 * t12 - t16 * t19) * t17);
    tzk0 = 0.2e1 * t24;

    # This function also can be used to calculate eps and Vxc

    t25 = params_beta * params_beta;
    t26 = 0.1e1 / t25;
    t27 = ρ * ρ;
    t28 = 0.1e1 / t27;
    t32 = my_piecewise3(t4, 0, -0.79577471545947667883e-1 * t16 * t26 * t28);
    tvrho0 = 0.2e1 * ρ * t32 + 0.2e1 * t24;
    #println("tvrho0 = ", tvrho0)

    return tvrho0

end