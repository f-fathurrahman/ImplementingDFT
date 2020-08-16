function KS_solve_TRDCM!(
    Ham::Hamiltonian, psi::Array{Float64,2};
    NiterMax=200, betamix=0.5,
    etot_conv_thr=1e-6,
    diag_func=diag_LOBPCG!,
)

    Npoints = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    dVol = Ham.grid.dVol

    Rhoe = zeros(Float64,Npoints)
    Rhoe_new = zeros(Float64,Npoints)
    
    calc_rhoe!( Ham, psi, Rhoe )
    println("integ Rhoe = ", sum(Rhoe)*dVol)
    update!( Ham, Rhoe )
    evals = zeros(Float64,Nstates)

    Y = zeros(Float64, Npoints, 3*Nstates)
    R = zeros(Float64, Npoints, Nstates)
    P = zeros(Float64, Npoints, Nstates)
    G = zeros(Float64, 3*Nstates, 3*Nstates)
    T = zeros(Float64, 3*Nstates, 3*Nstates)
    B = zeros(Float64, 3*Nstates, 3*Nstates)
    A = zeros(Float64, 3*Nstates, 3*Nstates)
    C = zeros(Float64, 3*Nstates, 3*Nstates)

    D = zeros(Float64,3*Nstates)

    set1 = 1:Nstates
    set2 = Nstates+1:2*Nstates
    set3 = 2*Nstates+1:3*Nstates
    set4 = Nstates+1:3*Nstates
    set5 = 1:2*Nstates

    MaxInnerSCF = 3
    MAXTRY = 10
    FUDGE = 1e-12
    SMALL = 1e-12

    ethr_evals_last = 1e-5
    ethr = 0.1

    Ham.energies.NN = calc_E_NN( Ham.atoms, Ham.pspots )

    Etot_old = 0.0

    Nconverges = 0

    for iter in 1:NiterMax
        #
        Hpsi = op_H(Ham, psi)
        Hsub = psi' * Hpsi * dVol  # Symmetrize Hsub?
        #
        # Calculate residual
        #
        R = Hpsi - psi * Hsub
        #
        # Precondition
        #
        for ist in 1:Nstates
            @views ldiv!( Ham.precKin, R[:,ist] )
        end
        #
        # Construct subspace (FIXME: work directly with Y?)
        #
        @views Y[:,set1] = psi[:,:]
        @views Y[:,set2] = R[:,:]
        if iter > 1
            @views Y[:,set3] = P[:,:]
        end
        #
        # Project kinetic and ionic potential
        #
        if iter > 1
            KY = -0.5*Ham.Laplacian*Y + op_V_Ps_loc(Ham, Y)
            T = Y' * KY * dVol
            B = Y' * Y * dVol
            # B = 0.5*( B + B' ) # Need this?
        else
            # only set5=1:2*Nstates is active for iter=1
            @views KY = -0.5*Ham.Laplacian * Y[:,set5] + op_V_Ps_loc( Ham, Y[:,set5] )
            @views T[set5,set5] = Y[:,set5]' * KY * dVol
            @views bb = Y[set5,set5]' * Y[set5,set5] * dVol
            @views B[set5,set5] = 0.5*( bb + bb' )
        end

        # Reset G
        if iter > 1
            G = Matrix(1.0I, 3*Nstates, 3*Nstates) #eye(3*Nstates)
        else
            G[set5,set5] = Matrix(1.0I, 2*Nstates, 2*Nstates)
        end

        @printf("TRDCM iter: %3d\n", iter)

        sigma = 0.0  # reset sigma to zero at the beginning of inner SCF iteration
        numtry = 0
        Etot0 = sum(Ham.energies)

        println("Etot0 = ", Etot0)

        for iterscf = 1:MaxInnerSCF
            
            #
            # Project Hartree, XC potential, and nonlocal pspot if any
            #
            V_loc = Ham.V_Hartree + Ham.V_XC
            
            if iter > 1
                yy = Y
            else
                yy = Y[:,set5]
            end

            if Ham.pspotNL.NbetaNL > 0
                VY = op_V_loc( V_loc, yy ) + op_V_Ps_nloc( Ham, yy )
            else
                VY = op_V_loc( V_loc, yy )
            end
                
            if iter > 1
                A = T + yy'*VY*dVol
                A = 0.5*(A + A')
            else
                aa = T[set5,set5] + yy'*VY*dVol
                A = 0.5*(aa + aa')
            end

            if iter > 1
                BG = B*G[:,1:Nstates]  # Nocc=Nstates?
                C = BG*BG'
                C = 0.5*(C + C')
            else
                BG = B[set5,set5]*G[set5,1:Nstates]
                cc = BG*BG'
                C[set5,set5] = 0.5*(cc + cc')
            end

            #
            # apply trust region if necessary
            #
            if abs(sigma) > SMALL # sigma is not zero
                println("Trust region is imposed")
                if iter > 1
                    D, G = eigen(A - sigma*C, B)
                else
                    D[set5], G[set5,set5] = eigen(A[set5,set5] - sigma*C[set5,set5], B[set5,set5])
                end
            else
                if iter > 1
                    D, G = eigen(A, B)
                else
                    D[set5], G[set5,set5] = eigen(A[set5,set5], B[set5,set5])
                end
            end
            #
            evals[:] = D[1:Nstates] .+ sigma
            #
            # update wavefunction
            #
            if iter > 1
                psi = Y*G[:,set1]
                ortho_sqrt!(psi, dVol)  # is this necessary ?
            else
                psi = Y[:,set5]*G[set5,set1]
                ortho_sqrt!(psi, dVol)
            end

            calc_rhoe!( Ham, psi, Rhoe )
            println("integ Rhoe = ", sum(Rhoe)*dVol)
            update!( Ham, Rhoe )

            # Calculate energies once again
            calc_energies!( Ham, psi )
            Etot = sum(Ham.energies)

            println("Etot = ", Etot)

            if Etot > Etot0

                @printf("TRDCM: %f > %f: Trust region will be imposed\n", Etot, Etot0)

                # Total energy is increased, impose trust region
                # Do this for all kspin

                if iter == 1
                    gaps = D[2:2*Nstates] - D[1:2*Nstates-1]
                    gapmax = maximum(gaps)
                else
                    gaps = D[2:3*Nstates] - D[1:3*Nstates-1]
                    gapmax = maximum(gaps)
                end
                gap0 = D[Nstates+1] - D[Nstates]

                while (gap0 < 0.9*gapmax) && (numtry < MAXTRY)
                    println("Increase sigma to fix gap0: numtry = ", numtry)
                    @printf("gap0 : %f < %f\n", gap0, 0.9*gapmax)
                    if abs(sigma) < SMALL # approx for sigma == 0.0
                        # initial value for sigma
                        sigma = 2*gapmax
                    else
                        sigma = 2*sigma
                    end
                    @printf("fix gap0: sigma = %18.10f\n", sigma)
                    #
                    if iter > 1
                        D, G = eigen(A - sigma*C, B)
                        gaps = D[2:2*Nstates] - D[1:2*Nstates-1]
                    else
                        D[set5], G[set5,set5] = eigen(A[set5,set5] - sigma*C[set5,set5], B[set5,set5])
                        gaps = D[2:3*Nstates] - D[1:3*Nstates-1]
                    end
                    gapmax = maximum(gaps)
                    gap0 = D[Nstates+1] - D[Nstates]
                    numtry = numtry + 1
                end

            end # if Etot > Etot0

            println("sigma = ", sigma)
            numtry = 0  # reset numtry for this while loop
            while (Etot > Etot0) &
                  #(abs(Etot-Etot0) > FUDGE*abs(Etot0)) &
                  (numtry < MAXTRY)
                @printf("Increase sigma part 2 try#%d: %f > %f ?\n", numtry, Etot, Etot0)
                #
                # update wavefunction
                #
                if iter > 1
                    psi = Y*G[:,set1]
                    ortho_sqrt!(psi, dVol)
                else
                    psi = Y[:,set5]*G[set5,set1]
                    ortho_sqrt!(psi, dVol)
                end

                calc_rhoe!( Ham, psi, Rhoe )
                update!( Ham, Rhoe )
                println("integ Rhoe = ", sum(Rhoe)*dVol)
            
                # Calculate energies once again
                calc_energies!( Ham, psi )
                Etot = sum(Ham.energies)
                #
                if Etot > Etot0
                    println("Increase sigma part 2")
                    if abs(sigma) > SMALL # sigma is not 0
                        sigma = 2*sigma
                    else
                        sigma = 1.2*gapmax
                    end
                    @printf("sigma = %f\n", sigma)
                    if iter > 1
                        D, G = eigen(A - sigma*C, B)
                    else
                        D[set5], G[set5,set5] = eigen(A[set5,set5] - sigma*C[set5,set5], B[set5,set5])
                    end
                end
                numtry = numtry + 1
            end # while

            Etot0 = Etot

        end # end of inner SCF

        # Calculate energies once again
        calc_energies!(Ham, psi)
        Etot = sum(Ham.energies)
        diffE = abs(Etot - Etot_old)
        @printf("TRDCM: %5d %18.10f %18.10e\n", iter, Etot, diffE)

        if diffE < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        if Nconverges >= 2
            @printf("\nTRDCM is converged in iter: %d\n", iter)
            break
        end

        Etot_old = Etot

        # No need to update potential, it is already updated in inner SCF loop
        if iter > 1
            P = Y[:,set4]*G[set4,set1]
        else
            P = Y[:,set2]*G[set2,set1]
        end

        flush(stdout)

    end

    @printf("\n")
    @printf("----------------------------\n")
    @printf("Final Kohn-Sham eigenvalues:\n")
    @printf("----------------------------\n")
    @printf("\n")
    for i in 1:Nstates
        @printf("%3d %18.10f\n", i, evals[i])
    end

    @printf("\n")
    @printf("----------------------------\n")
    @printf("Final Kohn-Sham energies:\n")
    @printf("----------------------------\n")
    @printf("\n")
    println(Ham.energies, banner=false)

    return
end