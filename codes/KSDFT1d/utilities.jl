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
