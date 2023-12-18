# https://discourse.julialang.org/t/how-to-generate-a-random-unitary-matrix-perfectly-in-julia/34102
function random_unitary_matrix(N::Int64)
    #x = (rand(N,N) + rand(N,N)*im) / sqrt(2)
    X = rand(Float64, N,N)
    F = qr(X)
    diagR = sign.(real(diag(F.R)))
    diagR[diagR.==0] .= 1
    diagRm = diagm(diagR)
    U = F.Q * diagRm
    return U
end

function transform_psi_Haux!(psi, Haux)
    λ, Urot = eigen(Hermitian(Haux))
    psi[:,:] = psi[:,:]*Urot
    Haux[:,:] = diagm(0 => λ)
    return Urot
end

function ortho_sqrt( psi::Array{Float64,2} )
    Udagger = inv(sqrt(psi'*psi))
    return psi*Udagger
end

function ortho_sqrt!( psi::Array{Float64,2} )
    Udagger = inv(sqrt(psi'*psi))
    psi[:,:] = psi*Udagger
    return
end

function generate_random_wavefunc(Ham)
    hx = Ham.grid.hx
    Npoints = Ham.grid.Npoints
    Nstates = Ham.electrons.Nstates
    psi = ortho_sqrt(rand(Npoints,Nstates))
    psi .*= (1.0/sqrt(hx))
    return psi
end

function print_matrix(A, Nrows, Ncols)
    for i in 1:Nrows
        for j in 1:Ncols
            @printf("%10.5f", A[i,j])
        end
        println()
    end
    return
end

function prec_invK(Ham::Hamiltonian1d, v)
    return inv(Ham.Kmat)*v
end


# Ideal preconditioner (expensive to calculate)
function prec_invHam(Ham::Hamiltonian1d, v)
    Kmat = Ham.Kmat
    Vtot = Ham.potentials.Total
    Hmat = Kmat + diagm( 0 => Vtot[:,1] )
    λ = eigvals(Hmat)
    vout = similar(v)
    for i in 1:size(v,2)
        @views vout[:,i] = inv(Hmat - λ[i]*I)*v[:,i]
    end
    return vout
end

