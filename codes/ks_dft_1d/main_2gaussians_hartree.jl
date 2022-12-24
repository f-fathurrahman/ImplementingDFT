using Printf
using LinearAlgebra

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("INC_sch_1d.jl")

# parameters?
function pot_gaussian( x, x0 )
    return -exp(-(x - x0)^2)/sqrt(1.0/π)^3
end

function calc_rhoe!(Focc, psi, rhoe)
    Npoints = size(psi, 1)
    Nstates = size(psi, 2)
    fill!(rhoe, 0.0)
    for ist in 1:Nstates
        for ip in 1:Npoints
            rhoe[ip] += Focc[ist] * psi[ip,ist] * psi[ip,ist]
        end
    end
    return
end


# a is in bohr, softening parameter of soft Coulomb potential
function Poisson_solve_sum!( xgrid, h,
    rho::Vector{Float64}, V::Vector{Float64}; a = 1.0
)
    Npoints = size(rho,1) 
    fill!(V, 0.0)
    for ip in 1:Npoints
        xi = xgrid[ip]
        for jp in 1:Npoints
            xj = xgrid[jp]
            dr = sqrt( (xi - xj)^2 + a^2 )
            V[ip] += rho[jp]/dr
        end
        V[ip] = V[ip]*h
    end
    return
end 

function calc_Ekin(D2, Focc, psi, h)
    Ekin = 0.0
    Nstates = size(psi, 2)
    Kpsi = -0.5*D2*psi
    for ist in 1:Nstates
        Ekin += Focc[ist]*dot(psi[:,ist], Kpsi[:,ist])*h
    end
    return Ekin
end


function main()
    # Initialize the grid points
    xmin = -5.0
    xmax =  5.0
    Npoints = 51
    xgrid, h = init_FD1d_grid(xmin, xmax, Npoints)
    
    # Build 2nd derivative matrix
    D2 = build_D2_matrix_9pt(Npoints, h)
    
    # Potential
    Vion = pot_gaussian.(xgrid, -1.0) + pot_gaussian.(xgrid, 1.0)
    
    #plt.clf()
    #plt.plot(xgrid, Vion)
    #plt.grid(true)
    #plt.savefig("IMG_Vion.png")

    Vtot = zeros(Float64, Npoints)
    Vhartree = zeros(Float64, Npoints)
    rhoe = zeros(Float64, Npoints)
    rhoe_new = zeros(Float64, Npoints)

    Nelectrons = 4
    Nstates = round(Int64, Nelectrons/2) # XXX only for even number of electrons
    Focc = 2*ones(Nstates)
    β_mix = 0.2

    Etot = Inf
    Etot_old = Etot
    Enn = 2*2/2 # Z_i*Z_j/r_ij

    for iter_scf in 1:100

        # Hamiltonian
        Vtot[:] = Vion[:] + Vhartree[:]
        Ham = -0.5*D2 + diagm( 0 => Vtot )
    
        # Solve the eigenproblem
        evals_all, evecs_all = eigen( Ham )    
        psi = evecs_all[:,1:Nstates]
        evals = evals_all[1:Nstates]
        # Renormalize
        psi[:] = psi[:]/sqrt(h)
        
        @printf("Eigenvalues\n")
        for ist in 1:Nstates
            @printf("%5d %18.10f\n", ist, evals[ist])
        end
        #println("Check normalization")
        #for ist in 1:Nstates
        #    println("dot(1,ist) = ", dot(psi[:,1], psi[:,ist])*h)
        #end

        calc_rhoe!(Focc, psi, rhoe_new)
        println("integ rhoe_new = ", sum(rhoe_new)*h)

        Ekin = calc_Ekin(D2, Focc, psi, h)
        Ehartree = 0.5*dot(rhoe_new, Vhartree)*h
        Eion = dot(rhoe_new, Vion)*h
        Etot =  Ekin + Ehartree + Eion + Enn

        ΔE = abs(Etot - Etot_old)
        mae_rhoe = sum(abs.(rhoe - rhoe_new))/Npoints
        @printf("%3d %18.10f %10.5e %10.5e\n", iter_scf, Etot, ΔE, mae_rhoe)

        if mae_rhoe < 1e-7
            println("Converged")
            Etot = sum(evals) + 0.5*dot(Vhartree, rhoe_new)*h + dot(rhoe_new,Vion)*h
            @printf("Ekin     = %18.10f\n", Ekin)
            @printf("Ehartree = %18.10f\n", Ehartree)
            @printf("Eion     = %18.10f\n", Eion)
            @printf("Ebands   = %18.10f\n", 2*sum(evals)) # FIXME factor of 2
            break
        end

        # Mix
        if iter_scf >= 2
            rhoe[:] = β_mix*rhoe_new[:] + (1 - β_mix)*rhoe[:]
        else
            rhoe[:] = rhoe_new[:]
        end
        Etot_old = Etot

        Poisson_solve_sum!(xgrid, h, rhoe, Vhartree)
        println("sum Vhartree = ", sum(Vhartree))

    end

end

main()