function calc_grad!(
    Ham, psi, g, Hsub
)
    ispin = 1 # FIXED

    Nstates = size(psi,2)
    Focc = Ham.electrons.Focc

    Hpsi = op_H( Ham, psi )
    Hsub[:,:] = psi' * Hpsi
    Hpsi = Hpsi - psi*Hsub

    for ist in 1:Nstates
        @views g[:,ist] = Focc[ist,ispin] * Hpsi[:,ist]
    end

    return
end

# Gradient for Haux
# The real input is actually stored in Ham.electrons.ebands which is calculated
# from  
function calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)

    Nspin = 1
    Nstates = Ham.electrons.Nstates

    fprime = zeros(Float64,Nstates)
    fprimeNum = zeros(Float64,Nstates)
    dmuNum = zeros(Float64,Nspin)
    dmuDen = zeros(Float64,Nspin)

    # These variables are not updated or calculated here
    # They are assumed to be calculated elsewhere
    ebands = Ham.electrons.ebands
    kT = Ham.electrons.kT
    E_fermi = Ham.electrons.E_fermi

    w = 2.0 # For Nspin=2
    ispin = 1 # FIXED

    for ist in 1:Nstates
        fprime[ist] = smear_fermi_prime( ebands[ist,ispin], E_fermi, kT )
        fprimeNum[ist] = fprime[ist] * ( real(Hsub[ist,ist]) - ebands[ist,ispin] )
    end
    dmuNum[ispin] = dmuNum[ispin] + w * sum(fprimeNum)
    dmuDen[ispin] = dmuDen[ispin] + w * sum(fprime)

    dmuContrib = sum(dmuNum)/sum(dmuDen)
    dBzContrib = 0.0 # not used
    println("fprime     = ", fprime)
    println("fprimeNum  = ", fprimeNum)
    println("dmuContrib = ", dmuContrib)

    gradF0 = zeros(Nstates,Nstates)
    gradF = zeros(Nstates,Nstates)

    g_tmp = zeros(Nstates,Nstates)

    gradF0[:,:] = Hsub - diagm( 0 => Ham.electrons.ebands[:,ispin] )
    gradF[:,:] = copy(gradF0)
    for ist in 1:Nstates
        gradF[ist,ist] = gradF0[ist,ist] - dmuContrib # FIXME: not tested for spinpol
    end
    g_tmp[:,:] = grad_smear( smear_fermi, smear_fermi_prime, ebands[:,ispin], E_fermi, kT, gradF )
    
    g_Haux[:,:] = w * 0.5 * (g_tmp' + g_tmp)
    Kg_Haux[:,:] = -copy(gradF0) #-0.1*copy(gradF0)

    return

end


# Eq (24)
function calc_dFdmu(Hsub, ebands, Focc_in, kT)
    Nstates = size(ebands,1)
    dFdmu = zeros(Nstates)
    Focc = 0.5*Focc_in # normalized to 1
    for ist in 1:Nstates
        dFdmu[ist] = (Hsub[ist,ist] - ebands[ist,1])*Focc[ist,1]*(1 - Focc[ist,1])
    end
    return dFdmu/kT
end


function offdiag_elements( Hsub, ebands, E_f::Float64, kT::Float64 )
    Nstates = size(evals, 1)
    mat = zeros(Nstates,Nstates)
    for j in 1:Nstates, i in 1:Nstates
        de = ebands[i] - ebands[j]
        if abs(de) > 1e-6
            mat[i,j] = Hsub[i,j] * ( smear_fermi(ebands[i], E_f, kT) - smear_fermi(ebands[j], E_f, kT) ) / de
        end
    end
    return mat
end