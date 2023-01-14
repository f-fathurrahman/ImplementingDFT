#function test01()

    Random.seed!(1234)

    Ham = init_Hamiltonian()

    hx = Ham.grid.hx
    Npoints = Ham.grid.Npoints
    Nelectrons = Ham.electrons.Nelectrons
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc

    psi = generate_random_wavefunc(Ham)

    Haux = psi' * (Ham*psi) * hx # Hsub, subspace Hamiltonian
    #Haux = rand(Nstates,Nstates)
    #Haux = 0.5*(Haux + Haux')

    update_from_wavefunc!(Ham, psi)
    update_from_Haux!(Ham, Haux)

    println(Ham.electrons.ebands)

    g = zeros(Npoints,Nstates)
    Hsub = zeros(Nstates,Nstates)

    calc_grad!(Ham, psi, g, Hsub)

    g_Haux = zeros(Nstates,Nstates)
    Kg_Haux = zeros(Nstates,Nstates)

    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)

#end


#=
Ham = init_Hamiltonian()
Ham.energies.NN = calc_E_NN(Ham.atoms)
psi2 = deserialize("TEMP_psi.dat")
evals2 = deserialize("TEMP_evals.dat")
Haux2 = diagm(0 => evals2[:,1])
=#
