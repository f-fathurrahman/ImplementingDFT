function PoissonYukawa1d_solve!(
    κ::Float64, ε0::Float64,
    gvec::GVectors1d,
    ρ_in::Vector{Float64},
    V_out::Vector{Float64}
)
    G2 = gvec.G2
    Ng = gvec.Ng # should be the same as Ns

    # Use Yukawa rather than the bare Coulomb potential.  The method does
    # not depend critically on the fact that the Coulomb kernel is to be
    # used.
    Vtemp = fft(ρ_in)
    Vtemp[1] = 0.0
    for ig in 2:Ng
        Vtemp[ig] *= 4π/(G2[ig] + κ^2) # multiply Vtemp
    end
    # Vtemp = 4*pi*invkx2.*fft(rho)  # we can decorate this one
    
    # invkx2 = [1/YukawaK^2;1./(kx(2:end).^2+YukawaK^2)];
    # The 1/YukawaK^2 factor is not important for neutral defect calculation, but
    # might be important for the charged defect calculation.

    V_out[:] = real( ifft(Vtemp) ) / ε0
    return
end