function calc_epsxc_1d( xc_calc::LibxcXCCalculator, Rhoe::Array{Float64,1} )

    Npoints = size(Rhoe)[1]
    Nspin = 1
    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    #Libxc_xc_func_init(ptr, 600, Nspin)  # LDA_X_1D_EXPONENTIAL
    Libxc_xc_func_init(ptr, 21, Nspin)  # LDA_X_1D_SOFT
    Libxc_xc_lda_exc!(ptr, Npoints, Rhoe, eps_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, 18, Nspin) # LDA_C_1D_CSC
    Libxc_xc_lda_exc!(ptr, Npoints, Rhoe, eps_c)
    Libxc_xc_func_end(ptr)

    #
    Libxc_xc_func_free(ptr)

    return eps_x + eps_c

end

function calc_Vxc_1d( xc_calc::LibxcXCCalculator, Rhoe::Array{Float64,1} )

    Npoints = size(Rhoe)[1]
    Nspin = 1
    v_x = zeros(Float64,Npoints)
    v_c = zeros(Float64,Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, 21, Nspin)
    Libxc_xc_lda_vxc!(ptr, Npoints, Rhoe, v_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, 18, Nspin)
    Libxc_xc_lda_vxc!(ptr, Npoints, Rhoe, v_c)
    Libxc_xc_func_end(ptr)

    #
    Libxc_xc_func_free(ptr)

    return v_x + v_c

end

