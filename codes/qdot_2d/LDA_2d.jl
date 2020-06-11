function calc_E_xc_V_xc_2d( rho::Array{Float64,1}; dVol=1.0 )
  
    Npoints = size(rho,1)
  
    Ex = 0.0
    Ec = 0.0
  
    Vx = zeros(Float64, Npoints)
    Vc = zeros(Float64, Npoints)
  
    for ip in 1:Npoints
        dex, Vx[ip] = calc_ex_vx_2d( rho[ip] )
        Ex = Ex + dex*rho[ip]
        dec, Vc[ip] = calc_ec_vc_2d( rho[ip] )
        Ec = Ec + dec*rho[ip]
    end

    Ex = Ex*dVol
    Ec = Ec*dVol
  
    E_xc = Ex + Ec
    V_xc = Vx + Vc

    return E_xc, V_xc
end


function calc_E_xc_2d( rho::Array{Float64,1}; dVol=1.0 )

    Npoints = size(rho,1)

    Ex = 0.0
    Ec = 0.0

    for ip in 1:Npoints
        dex = calc_ex_2d( rho[ip] )
        Ex = Ex + dex*rho[ip]
        dec = calc_ec_2d( rho[ip] )
        Ec = Ec + dec*rho[ip]
    end

    Ex = Ex*dVol
    Ec = Ec*dVol
    E_xc = Ex + Ec

    return E_xc
end

function calc_V_xc_2d( rho::Array{Float64,1}; dVol=1.0 )
    Npoints = size(rho,1)

    Vx = zeros(Float64, Npoints)
    Vc = zeros(Float64, Npoints)
    
    for ip in 1:Npoints
        Vx[ip] = calc_vx_2d( rho[ip] )
        Vc[ip] = calc_vc_2d( rho[ip] )
    end
    
    V_xc = Vx + Vc
    
    return V_xc
end


function calc_ex_vx_2d( ρ_ )
    a_x = -1.06384608107049
    SMALL = 1.0e-15
    ρ = max(0.0, ρ_)
    if ρ < SMALL
        return 0.0, 0.0
    end
    ex = a_x*sqrt(ρ)
    vx = 1.5*ex
    return ex, vx
end

function calc_ex_2d( ρ_ )
    a_x = -1.06384608107049
    SMALL = 1.0e-15
    ρ = max(0.0, ρ_)
    if ρ < SMALL
        return 0.0
    end
    ex = a_x*sqrt(ρ)
    return ex
end

function calc_vx_2d( ρ_ )
  a_x = -1.06384608107049
  SMALL = 1.0e-15
  ρ = max(0.0, ρ_)
  if ρ < SMALL
      return 0.0
  end
  ex = a_x*sqrt(ρ)
  vx = 1.5*ex
  return vx
end


function calc_ec_vc_2d( ρ_ )
  
    a = -0.1925
    b = 0.0863136
    c = 0.0572384
    e = 1.0022
    f = -0.02069
    g = 0.33997
    h = 1.747e-2
    d = -a*h

    SMALL = 1.0e-15
    ρ = max(0.0, ρ_)
    if ρ < SMALL
      return 0.0, 0.0
    end
    rs = sqrt(1.0/(pi*ρ))

    ec = a + (b*rs + c*rs^2 + d*rs^3) * log(1.0 + 1.0/(e*rs + f*sqrt(rs)^3 + g*rs^2 + h*rs^3) )

    efe  = e*rs + f*sqrt(rs)^3 + g*rs^2 + h*rs^3
    efep = e + 1.5*f*sqrt(rs) + 2.0*g*rs + 3.0*h*rs^2
    lg = log(1.0 + 1.0/efe)
    x  = ( (b*rs + c*rs^2 + d*rs^3)*efep )/(efe^2 + efe)
    dalphadrs = (b + 2.0*c*rs + 3.0*d*rs^2)*lg - x
    vc = ec - 0.5*rs*dalphadrs

    return ec, vc
end


function calc_ec_2d( ρ_ )
  
    a = -0.1925
    b = 0.0863136
    c = 0.0572384
    e = 1.0022
    f = -0.02069
    g = 0.33997
    h = 1.747e-2
    d = -a*h

    SMALL = 1.0e-15
    ρ = max(0.0, ρ_)
    if ρ < SMALL
      return 0.0
    end
    rs = sqrt(1.0/(pi*ρ))

    ec = a + (b*rs + c*rs^2 + d*rs^3) * log(1.0 + 1.0/(e*rs + f*sqrt(rs)^3 + g*rs^2 + h*rs^3) )

    return ec
end


function calc_vc_2d( ρ_ )
  
    a = -0.1925
    b = 0.0863136
    c = 0.0572384
    e = 1.0022
    f = -0.02069
    g = 0.33997
    h = 1.747e-2
    d = -a*h

    SMALL = 1.0e-15
    ρ = max(0.0, ρ_)
    if ρ < SMALL
      return 0.0
    end
    rs = sqrt(1.0/(pi*ρ))

    ec = a + (b*rs + c*rs^2 + d*rs^3) * log(1.0 + 1.0/(e*rs + f*sqrt(rs)^3 + g*rs^2 + h*rs^3) )

    efe  = e*rs + f*sqrt(rs)^3 + g*rs^2 + h*rs^3
    efep = e + 1.5*f*sqrt(rs) + 2.0*g*rs + 3.0*h*rs^2
    lg = log(1.0 + 1.0/efe)
    x  = ( (b*rs + c*rs^2 + d*rs^3)*efep )/(efe^2 + efe)
    dalphadrs = (b + 2.0*c*rs + 3.0*d*rs^2)*lg - x
    vc = ec - 0.5*rs*dalphadrs

    return vc
end


#=
# Correlation term missing
function fxc_lda(n, fxc)

  real(8), parameter :: a_x = -1.06384608107049d0
  SMALL = 1.0d-15

  ! Sanity check
  dens = max(0.0d0, n)

  ! If the density is too small, return zero  
  if(dens < SMALL) then
    fxc = 0.0d0
    return
  endif

  fxc = 0.75d0*a_x/sqrt(dens)

end
=#

