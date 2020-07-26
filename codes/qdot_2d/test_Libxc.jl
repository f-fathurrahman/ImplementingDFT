using Libxc

function main()
    x_func = Functional(:lda_x_2d)
    c_func = Functional(:lda_c_2d_amgb)

    rho = [1.0, 2.0, 3.0, 4.0]
    epsxc = evaluate(x_func, rho=rho, derivatives=0).zk +
            evaluate(c_func, rho=rho, derivatives=0).zk

    println("epsxc = ", epsxc)
    println("Exc = ", sum(rho .* epsxc))
end

main()