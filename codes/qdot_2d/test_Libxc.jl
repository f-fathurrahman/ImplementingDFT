using Libxc

function test_LDA()
    x_func = Functional(:lda_x_2d)
    c_func = Functional(:lda_c_2d_amgb)
    #
    rho = [1.0, 2.0, 3.0, 4.0]
    epsxc = evaluate(x_func, rho=rho, derivatives=0).zk +
            evaluate(c_func, rho=rho, derivatives=0).zk
    #
    Vxc = evaluate(x_func, rho=rho, derivatives=1).vrho +
          evaluate(c_func, rho=rho, derivatives=1).vrho
    #
    println("epsxc = ", epsxc)
    println("Exc = ", sum(rho .* epsxc))
    println("Vxc = ", Vxc)
end
#test_LDA()

# Using inplace evaluate!
function test_LDA_spin()
    x_func = Functional(:lda_x_2d, n_spin=2)
    c_func = Functional(:lda_c_2d_amgb, n_spin=2)
    #
    rho = [1.0, 2.0, 3.0, 4.0]
    rho2 = zeros(Float64, 2*length(rho))
    rho2[1:2:end] = rho2[2:2:end] = rho ./ 2
    #
    eps_x = zeros(Float64, length(rho))
    evaluate!(x_func, rho=rho2, zk=eps_x)
    println("eps_x = ", eps_x)
    #
    eps_c = zeros(Float64, length(rho))
    evaluate!(c_func, rho=rho2, zk=eps_c)
    println("eps_c = ", eps_c)

    #evaluate!(c_func, rho=rho2, derivatives=0).zk
    #
    #
    #Vxc = zeros(length(Vxc))
    #evaluate!(x_func, rho=rho2, derivatives=1, vrho=Vxc)
    #evaluate!(c_func, rho=rho2, derivatives=1).vrho
    #
    #println("Exc = ", sum(rho .* epsxc))
    #println("Vxc = ", Vxc)
end
test_LDA_spin()