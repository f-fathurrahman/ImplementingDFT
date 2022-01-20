function hartree_pot_bc(rho, H::Hamiltonian)
    # computes the Hartree potential and energy in periodic boundary potential
    # by solving the 1-d Yakuwa equation.
    # we call the (hopefully compiler optmized version)
    return hartree_pot_bc(rho, H.Ls, H.YukawaK, H.epsil0)
end

function hartree_pot_bc(rho::Array{Float64,2}, Ls::Float64, YukawaK::Float64, epsil0::Float64)
    # Calculate the Hartree potential and the energy in periodic boundary potential
    # this function is called from the one in the Ham. This is anotated to give the jit
    # compiler more information
    Ns_glb = length(rho);
    kx = 2*pi* vcat(collect(0:Ns_glb/2), collect(-Ns_glb/2+1:-1))./Ls
    # Use Yukawa rather than the bare Coulomb potential.  The method does
    # not depend critically on the fact that the Coulomb kernel is to be
    # used.
    invkx2 = [0;1 ./(kx[2:end].^2 .+ YukawaK^2)]
    # invkx2 = [1/YukawaK^2;1./(kx(2:end).^2+YukawaK^2)];
    # The 1/YukawaK^2 factor is not important for neutral defect calculation, but
    # might be important for the charged defect calculation.
    Vtemp = 4*pi*invkx2.*fft(rho)  # we can decorate this one
    V = real(ifft(Vtemp))
    return V / epsil0
end



############################################################################
##  Trying different optimizations

function hartree_pot_bc_opt(rho::Array{Float64,2}, Ls::Float64,
                            YukawaK::Float64, epsil0::Float64)
    # Calculate the Hartree potential and the energy in periodic boundary potential
    # this function is called from the one in the Ham. This is anotated to give the jit
    # compiler more information
    rhoFourier = rfft(rho)
    inv_yukawa_fourier_mult!(rhoFourier,Ls,YukawaK)
    V  = irfft(rhoFourier, size(rho,1));
    return V / epsil0;
end


function inv_yukawa_fourier_mult!(R::Array{ComplexF64,2}, Ls::Float64,
                                          YukawaK::Float64 )
    c = (2 * pi / Ls)^2
    c1 = 4*pi
    R[1] = 0
    @inbounds @simd for ii = 2:size(R)[1]
        R[ii] = c1*R[ii]/((ii-1)^2*c + YukawaK^2);
    end
end


# using 1D vector in general should produce faster resutls

function hartree_pot_bc_opt_vec(rho::Vector{Float64}, Ls::Float64,
                                YukawaK::Float64, epsil0::Float64)
    # Calculate the Hartree potential and the energy in periodic boundary potential
    # this function is called from the one in the Ham. This is anotated to give the jit
    # compiler more information
    rhoFourier = rfft(rho)
    inv_yukawa_fourier_mult_vec!(rhoFourier,Ls,YukawaK)
    V  = irfft(rhoFourier, size(rho,1))/ epsil0;
    return V;
end


function inv_yukawa_fourier_mult_vec!(
    R::Vector{ComplexF64}, Ls::Float64,
    YukawaK::Float64
)
    c = (2 * pi / Ls)^2
    c1 = 4*pi
    R[1] = 0
    @inbounds @simd for ii = 2:size(R)[1]
        R[ii] = c1*R[ii]/((ii-1)^2*c + YukawaK^2);
    end
end
