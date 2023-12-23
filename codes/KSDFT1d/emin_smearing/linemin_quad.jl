
function linemin_quad_psi(Ham, psi, Haux, g, g_Haux, d, d_Haux, E1)
    hx = Ham.grid.hx
    #
    α_t = 1e-5
    psic = psi + α_t*d
    Hauxc = deepcopy(Haux)
    #
    prepare_psi_Haux!(psic, Hauxc, hx)
    Etrial = calc_Lfunc_Haux!(Ham, psic, Hauxc)
    ΔEdir = 2*dot(g, α_t*d)*hx
    #
    println("linemin_quad_psi: E1     = ", E1)
    println("linemin_quad_psi: Etrial = ", Etrial)
    if Etrial > E1
        println("linemin_quad_psi: WARNING!!! Etrial is higher than E1")
        return α_t, false
    end
    println("linemin_quad_psi: ΔEdir  = ", ΔEdir)
    println("linemin_quad_psi: ratio of energy diff = ", (Etrial - E1)/ΔEdir)
    #
    c = ( Etrial - (E1 + ΔEdir) ) / α_t
    println("linemin_quad_psi: c = ", c)
    return abs(ΔEdir/(2*c)), true
end


function linemin_quad_Haux(Ham, psi, Haux, g, g_Haux, d, d_Haux, E1)
    hx = Ham.grid.hx
    #
    α_t = 1e-5
    psic = deepcopy(psi)
    Hauxc = Haux + α_t*d_Haux
    #
    prepare_psi_Haux!(psic, Hauxc, hx)
    Etrial = calc_Lfunc_Haux!(Ham, psic, Hauxc)
    ΔEdir = dot(g_Haux, α_t*d_Haux)
    #
    println("linemin_quad_Haux: E1     = ", E1)
    println("linemin_quad_Haux: Etrial = ", Etrial)
    if Etrial > E1
        println("linemin_quad_Haux: WARNING!!! Etrial is higher than E1")
        return α_t, false
    end
    println("linemin_quad_Haux: ΔEdir  = ", ΔEdir)
    println("linemin_quad_Haux: ratio of energy diff = ", (Etrial - E1)/ΔEdir)
    #
    c = ( Etrial - (E1 + ΔEdir) ) / α_t
    println("linemin_quad_Haux: c = ", c)
    return abs(ΔEdir/(2*c)), true
end



# Vary psi and Haux simulataneously
function linemin_quad(Ham, psi, Haux, g, g_Haux, d, d_Haux, E1)
    hx = Ham.grid.hx
    #
    α_t = 1e-4
    psic = psi + α_t*d # trial wavefunc
    Hauxc = Haux + α_t*d_Haux
    #
    Udagger = inv(sqrt(psic'*psic)) ./ sqrt(hx) # rotation
    psic[:,:] = psic*Udagger # orthonormalize
    Urot2 = transform_psi_Haux!(psic, Hauxc) # make Haux diagonal 
    #
    Etrial = calc_Lfunc_Haux!(Ham, psic, Hauxc)
    ΔEdir = 2*dot(g, α_t*d)*hx + dot(g_Haux, α_t*d_Haux)
    #
    println("LineMin: E1     = ", E1)
    println("LineMin: Etrial = ", Etrial)
    if Etrial > E1
        println("LineMin: WARNING!!! Etrial is higher than E1")
        return α_t, false
    end
    println("LineMin: ΔEdir  = ", ΔEdir)
    println("LineMin: ratio of energy diff = ", (Etrial - E1)/ΔEdir)
    #
    c = ( Etrial - (E1 + ΔEdir) ) / α_t
    println("LineMin: c = ", c)
    #if c < 0.0
    #    println("LineMin: Negative curvature, returning α_t")
    #    return α_t, true
    #end
    return abs(ΔEdir/(2*c)), true
end