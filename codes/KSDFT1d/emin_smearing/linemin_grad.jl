function linemin_grad(Ham, psi, Haux, g, g_Haux, d, d_Haux, E1)
    #
    hx = Ham.grid.hx
    #
    gt = zeros(Float64, size(psi))
    gt_Haux = zeros(Float64, size(Haux))
    Kgt_Haux = zeros(Float64, size(Haux))
    Hsub = zeros(Float64, size(Haux))
    #
    α_t = 1e-3
    psic = psi + α_t*d # trial wavefunc
    Hauxc = Haux + α_t*d_Haux
    #
    Udagger = inv(sqrt(psic'*psic)) ./ sqrt(hx) # rotation
    psic[:,:] = psic*Udagger # orthogonalize
    Hauxc = Udagger' * Hauxc * Udagger # rotate Haux
    Urot2 = transform_psi_Haux!(psic, Hauxc) # make Haux diagonal 
    #
    Etrial = calc_Lfunc_Haux!(Ham, psic, Hauxc) # need to call this?
    calc_grad_Lfunc_Haux!(Ham, psic, Hauxc, gt, Hsub, gt_Haux, Kgt_Haux)
    println("LineMin: Etrial = ", Etrial)
    if Etrial > E1
        println("LineMin: Etrial is not smaller")
    end
    #
    denum1 = 2*sum(conj(g-gt).*d)*hx + sum(conj(g_Haux-gt_Haux).*d_Haux)
    println("LineMin: denum1 = ", denum1)
    if denum1 != 0.0 # check for small
        num1 = 2*sum(conj(g).*d)*hx + sum(conj(g_Haux).*d_Haux)
        println("num1 = ", num1)
        α = α_t*num1/denum1
    else
        α = 0.0
    end
    println("LineMin: α = ", α)
    return α
end