function test_hartree_pot_bc(Ham)
    #
    Ns = Ham.Ns
    rho_in = zeros(Float64, Ns, 1)
    fill!( rho_in, 1.0 )
    rho_in[2] = 2.0
    VHartree = hartree_pot_bc(rho_in, Ham)

    println("Test hartree_pot_bc, Ns = ", Ns)
    for ip in 1:10
        @printf("%3d %18.10f\n", ip, VHartree[ip])
    end
end