function XC_c_vwn( Rhoe )

    third = 1.0/3.0
    pi34 = 0.6203504908994  # pi34=(3/4pi)^(1/3)
    rs = pi34/Rhoe^third
    rs12 = sqrt(rs)

    A = 0.0310907 # In Hartree
    b = 3.72744
    c = 12.9352
    x0 = -0.10498

    Q = sqrt(4.0 * c - b * b)
    f1 = 2.0 * b / Q
    f2 = b * x0 / (x0 * x0 + b * x0 + c)
    f3 = 2.0 * (2.0 * x0 + b) / Q
    fx = rs + b * rs12 + c
    qx = atan( Q/(2.0 * rs12 + b) )
    # The correlation energy density
    ec = A * ( log(rs/fx) + f1 * qx - f2 * ( log( (rs12 - x0)^2 /fx) + f3 * qx) )
    
    tx = 2.0 * rs12 + b
    tt = tx * tx + Q * Q
    # The potential
    vc = ec - rs12 * A / 6.0 * (2.0 / rs12 - tx / fx - 4.0 * b / tt -
         f2 * (2.0 / (rs12 - x0) - tx / fx - 4.0 * (2.0 * x0 + b) / tt) )

    return ec, vc

end