# Create the pseudo-charge for Coulomb interaction
# $\rho_a = \sum_{N_{\text{atoms}}} - Z_j/\sqrt{2 \pi \sigma_i}$
function pseudocharge(gridpos, Ls, atoms, YukawaK, epsil0)
    
    Nsglb = length(gridpos)
    Natoms = atoms.Natoms

    Z     = atoms.Z
    sigma = atoms.sigma
    R     = atoms.R

    # Calculate the total pseudo-charge
    # NOTE: The pseudo-charge should not extend over twice the domain size!
    rhoa = zeros(Nsglb)
    for ia in 1:Natoms
        d = R[ia] .- gridpos
        d = d - round.(Integer, d/Ls)*Ls
        # Note: Z has minus sign in front of it!
        rhoa = rhoa - Z[ia]./sqrt(2*pi*sigma[ia]^2) .* (exp.(-0.5*(d./sigma[ia]).^2 ))
    end

    drhoa = zeros(Nsglb, Natoms)
    for j in 1:Natoms
        d = R[j] .- gridpos
        d = d - round.(Integer, d/Ls)*Ls
        drhoa[:,j1] = drhoa[:,j] - Z[j] ./sqrt(2*pi*sigma[j]^2)./sigma[j]^2 .* (exp.(-0.5*(d./sigma[j]).^2)).*d
    end
    return rhoa, drhoa
end
