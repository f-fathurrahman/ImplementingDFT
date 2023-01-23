mutable struct ElecVars
    psi::Matrix{Float64} # tall matrix: (Nbasis,Nstates)
    Hsub::Matrix{Float64} # square matrix: (Nstates,Nstates)
    Hsub_eigs::Vector{Float64} # column vector: (Nstates)
end

function ElecVars( Ham::Hamiltonian1d )
    return ElecVars( Ham, generate_random_wavefunc(Ham) )
end

function ElecVars( Ham::Hamiltonian1d, psi::BlochWavefunc )
    
    Nstates = Ham.electrons.Nstates
    hx = Ham.grid.hx
    Nspin = 1 # FIXED


    Hsub = zeros(Float64,Nstates,Nstates)
    Hsub_eigs = zeros(Float64,Nstates)

    update_from_wavefunc!(Ham, psi)    

    Hsub = psi' * (Ham * psi) * hx
    Hsub_eigs = eigvals(Hermitian(Hsub[i]))  # set Haux_eigs to eigenvalues of Hsub

    return ElecVars(psiks, Hsub, Hsub_eigs)
end