using OffsetArrays

function calc_eps_c_1d_csc(ρ)

    # exponentially screened interaction
    #params_orig = [  4.66,  0.0,  2.092, 3.735, 0.0, 1.379, 2.0, 23.63,  109.9,    1.837 ];
    # soft-Coulomb interaction
    #params_orig = [ 7.40, 1.120, 1.890, 0.0964,  0.0250,   2.0, 3.0, 2.431, 0.0142, 2.922];
    #params_orig = [ 32.1,   0.0,  3.77,  7.576, 0.0, 0.941, 2.0,  1.63,    0.198,  4.086]; # interaction=0, para idx = 4 (bb==1.0)

    params_orig = [18.40, 0.0,   7.501, 0.10185, 0.012827, 2.0, 3.0, 1.511, 0.258,  4.424] # bb=1.0 interaction=1

    params_para = OffsetArray(params_orig, 0:9)

    t1 = 0.1e1 / ρ;

    t2 = t1 / 0.2e1;

    t3 = params_para[4];

    t4 = ρ * ρ;

    t5 = 0.1e1 / t4;
    t8 = t2 + t3 * t5 / 0.4e1;

    t9 = params_para[7];
    t13 = params_para[9];
    t14 = t2^t13;
    t15 = params_para[8] * t14;
    t16 = 0.1e1 + t9 * t1 / 0.2e1 + t15;
    t17 = log(t16);
    t18 = t8 * t17;
    t21 = params_para[1];
    t24 = params_para[5];
    t25 = t2^t24;
    t26 = params_para[2] * t25;
    t29 = params_para[6];
    t30 = t2^t29;
    t31 = params_para[3] * t30;
    t33 = t21 * t1 + 0.2e1 * t26 + 0.2e1 * t31 + 0.2e1 * params_para[0];
    t34 = 0.1e1 / t33;
    tzk0 = -t18 * t34;

    return tzk0

end


function calc_V_c_1d_csc(ρ::Float64)

    params_orig = [18.40, 0.0,   7.501, 0.10185, 0.012827, 2.0, 3.0, 1.511, 0.258,  4.424] # bb=1.0 interaction=1

    params_para = OffsetArray(params_orig, 0:9)

    t1 = 0.1e1 / ρ;
    t2 = t1 / 0.2e1;
    t3 = params_para[4];
    t4 = ρ * ρ;
    t5 = 0.1e1 / t4;
    t8 = t2 + t3 * t5 / 0.4e1;
    t9 = params_para[7];
    t13 = params_para[9];
    t14 = t2^t13;
    t15 = params_para[8] * t14;
    t16 = 0.1e1 + t9 * t1 / 0.2e1 + t15;
    t17 = log(t16);
    t18 = t8 * t17;
    t21 = params_para[1];
    t24 = params_para[5];
    t25 = t2^t24;
    t26 = params_para[2] * t25;
    t29 = params_para[6];
    t30 = t2^t29;
    t31 = params_para[3] * t30;
    t33 = t21 * t1 + 0.2e1 * t26 + 0.2e1 * t31 + 0.2e1 * params_para[0];
    t34 = 0.1e1 / t33;
    tzk0 = -t18 * t34;

    # tzk0 is eps

    t37 = 0.1e1 / t4 / ρ;
    t40 = -t3 * t37 / 0.2e1 - t5 / 0.2e1;
    t41 = ρ * t40;
    t42 = t17 * t34;
    t44 = ρ * t8;
    t49 = -t9 * t5 / 0.2e1 - t15 * t13 * t1;
    t50 = 0.1e1 / t16;
    t52 = t49 * t50 * t34;
    t54 = t33 * t33;
    t55 = 0.1e1 / t54;
    t56 = t17 * t55;
    t64 = -0.2e1 * t26 * t24 * t1 - 0.2e1 * t31 * t29 * t1 - t21 * t5;
    t65 = t56 * t64;
    tvrho0 = -t41 * t42 - t44 * t52 + t44 * t65 + tzk0;

    return tvrho0
end