using LinearAlgebra, FFTW, Statistics, Random
using Serialization

# ============================================================================
# 1D grid and helper functions
# ============================================================================

struct Grid1D
    x::Vector{Float64}
    dx::Float64
    L::Float64
    N::Int
end

function Grid1D(; L=20.0, N=100)
    # N points, endpoints -L/2 and L/2 are both included
    dx = L / (N - 1)
    x = range(-L/2, L/2, length=N) |> collect
    Grid1D(x, dx, L, N)
end

dot_grid(x::Vector{Float64}, y::Vector{Float64}, dx::Float64) = sum(x .* y) * dx
dot_grid(x::AbstractMatrix, y::AbstractMatrix, dx::Float64) = (x' * y) * dx

function orthonormalize!(psi::Matrix{Float64}, dx::Float64)
    Nstates = size(psi, 2)
    for j in 1:Nstates
        for i in 1:j-1
            proj = dot_grid(psi[:,i], psi[:,j], dx)
            psi[:,j] .-= proj * psi[:,i]
        end
        norm2 = dot_grid(psi[:,j], psi[:,j], dx)
        norm2 < 1e-12 && error("norm too small in orthonormalization")
        psi[:,j] ./= sqrt(norm2)
    end
    return psi
end

# ============================================================================
# Hamiltonian operations (finite‑difference)
# ============================================================================

function apply_kinetic(psi::Matrix{Float64}, grid::Grid1D)
    N, Nstates = size(psi)
    dx = grid.dx
    invdx2 = 1.0 / dx^2
    Tpsi = similar(psi)
    for i in 1:N
        if i == 1
            # left boundary: ψ[0] = 0
            Tpsi[i,:] = -0.5 * invdx2 * (-2psi[i,:] + psi[i+1,:])
        elseif i == N
            # right boundary: ψ[N+1] = 0
            Tpsi[i,:] = -0.5 * invdx2 * (psi[i-1,:] - 2psi[i,:])
        else
            Tpsi[i,:] = -0.5 * invdx2 * (psi[i-1,:] - 2psi[i,:] + psi[i+1,:])
        end
    end
    return Tpsi
end



function apply_H(psi::Matrix{Float64}, V_tot::Vector{Float64}, grid::Grid1D)
    return apply_kinetic(psi, grid) .+ V_tot .* psi
end

# ============================================================================
# External potential (soft‑Coulomb)
# ============================================================================

function external_potential(grid::Grid1D, atoms::Vector{Float64}, Z::Vector{Float64}, a_soft::Float64)
    V = zeros(grid.N)
    for (pos, charge) in zip(atoms, Z)
        V .+= -charge ./ sqrt.((grid.x .- pos).^2 .+ a_soft^2)
    end
    return V
end

# ============================================================================
# Hartree potential via direct integration (O(N²))
# ============================================================================

function hartree_potential(ρ::Vector{Float64}, grid::Grid1D, a_soft::Float64)
    N = grid.N
    dx = grid.dx
    x = grid.x
    V_H = zeros(N)
    for i in 1:N
        s = 0.0
        for j in 1:N
            d = x[i] - x[j]
            s += ρ[j] / sqrt(d^2 + a_soft^2) * dx
        end
        V_H[i] = s
    end
    return V_H
end

# ============================================================================
# Fermi‑Dirac occupation and entropy (with degeneracy g)
# ============================================================================

fermi(ε::Float64, μ::Float64, kT::Float64, g) = g / (1 + exp((ε - μ) / kT))

function find_mu(ε::Vector{Float64}, N_electrons::Float64, kT::Float64, g; tol=1e-12)
    μ_min = minimum(ε) - 10kT
    μ_max = maximum(ε) + 10kT
    f_sum(μ) = sum(fermi.(ε, μ, kT, g)) - N_electrons
    f_min = f_sum(μ_min)
    f_max = f_sum(μ_max)
    f_min * f_max > 0 && error("Cannot bracket chemical potential")
    for _ in 1:200
        μ = (μ_min + μ_max) / 2
        f = f_sum(μ)
        abs(f) < tol && return μ
        if f * f_min > 0
            μ_min = μ
            f_min = f
        else
            μ_max = μ
            f_max = f
        end
    end
    error("find_mu did not converge")
end

function entropy(f::Vector{Float64}, g)
    S = 0.0
    for fi in f
        if fi > 0 && fi < g
            # -k [ f ln(f/g) + (g-f) ln(1 - f/g) ]
            p = fi / g
            S -= fi * log(p) + (g - fi) * log(1 - p)
        end
    end
    return S   # note: k is not included here; it will be multiplied later
end

# ============================================================================
# Free energy and related quantities
# ============================================================================

function calc_density(psi::Matrix{Float64}, f::Vector{Float64})
    N, Nstates = size(psi)
    ρ = zeros(N)
    for n in 1:Nstates
        ρ .+= f[n] * psi[:,n].^2
    end
    return ρ
end

function calc_free_energy(psi, f, Hpsi, Hmat, density, V_H, grid, V_ext, kT, g)
    dx = grid.dx
    Nstates = size(f, 1)
    E_kin = 0.0
    for ist in 1:Nstates
        E_kin += f[n] * Hmat[n,n]
    end
    E_ext = sum(density .* V_ext) * dx
    E_H = 0.5 * sum(density .* V_H) * dx
    S = entropy(f, g)
    println("E_kin  = $(E_kin)")
    println("E_ext  = $(E_ext)")
    println("E_kext = $(E_kin+E_ext)")
    println("E_H    = $(E_H)")
    println("mTS    = $(-kT*S)")
    return E_kin + E_ext + E_H - kT * S
end

# ============================================================================
# Gradients (adapted for degeneracy g)
# ============================================================================

function gradient_psi(psi, Hpsi, f, dx)
    Hmat = (psi' * Hpsi) * dx
    g_psi = similar(psi)
    for n in 1:size(psi,2)
        proj = psi * Hmat[:,n]
        g_psi[:,n] = f[n] * (Hpsi[:,n] - proj)
    end
    return g_psi, Hmat
end

function gradient_eta(psi, Hmat, f, ε, kT, g)
    Nstates = length(ε)
    g_eta = zeros(Nstates, Nstates)
    w = f .* (1.0 .- f / g)          # derivative factor: f(1 - f/g)
    S = sum(w)
    avg = S > 0 ? sum((diag(Hmat) - ε) .* w) / S : 0.0
    for n in 1:Nstates
        g_eta[n,n] = -(1/kT) * w[n] * (Hmat[n,n] - ε[n] - avg)
        for m in n+1:Nstates
            Δε = ε[n] - ε[m]
            if abs(Δε) > 1e-8
                g_eta[n,m] = (f[n] - f[m]) / Δε * Hmat[n,m]
                g_eta[m,n] = g_eta[n,m]
            end
        end
    end
    return g_eta
end

# ============================================================================
# Preconditioner for wavefunctions (FFT‑based, periodic)
# ============================================================================

function preconditioner_psi(grid::Grid1D; ε_cut=0.5)
    N = grid.N
    L = grid.L
    freqs = fftfreq(N, 2π/L)
    Tk = 0.5 * freqs.^2
    return 1.0 ./ max.(Tk, ε_cut)
end

function apply_preconditioner(r::Vector{Float64}, precond::Vector{Float64})
    r_fft = fft(r)
    r_fft .*= precond
    return real(ifft(r_fft))
end

function apply_preconditioner(R::Matrix{Float64}, precond::Vector{Float64})
    N, Nstates = size(R)
    res = similar(R)
    for n in 1:Nstates
        res[:,n] = apply_preconditioner(R[:,n], precond)
    end
    return res
end

# ============================================================================
# Main conjugate‑gradient minimization (with degeneracy)
# ============================================================================

function minimize_dft(;
    L=20.0, Ngrid=100,
    atoms=[0.0], Z=[1.0],
    a_soft_ee=1.0, a_soft_en=1.0,
    N_electrons=1.0,
    kT=0.01,
    Nstates=4,
    degeneracy=1,          # spin degeneracy (1 = spin‑polarised, 2 = unpolarised)
    κ=0.1,
    maxiter=200, tol=1e-10,
    verbose=true,
    seed=1234,
    psi=nothing,
    ε=nothing
)
    if seed !== nothing
        Random.seed!(seed)
    end

    grid = Grid1D(L=L, N=Ngrid)
    V_ext = external_potential(grid, atoms, Z, a_soft_en)
    precond_psi = preconditioner_psi(grid)

    # initial wavefunctions
    if isnothing(psi)
        psi = randn(Ngrid, Nstates)
    end
    orthonormalize!(psi, grid.dx)

    # initial diagonal η (guess energies)
    if isnothing(ε)
        ε = collect(range(-1.0, 1.0, length=Nstates))
    end
    μ = find_mu(ε, N_electrons, kT, degeneracy)
    f = fermi.(ε, μ, kT, degeneracy)

    ρ = calc_density(psi, f)
    V_H = hartree_potential(ρ, grid, a_soft_ee)
    V_tot = V_ext + V_H
    Hpsi = apply_H(psi, V_tot, grid)
    Hmat = dot_grid(psi, Hpsi, grid.dx)
    F = calc_free_energy(psi, f, Hpsi, Hmat, ρ, V_H, grid, V_ext, kT, degeneracy)

    X_psi = zeros(Ngrid, Nstates)
    X_eta = zeros(Nstates, Nstates)
    dot_old = 0.0
    eta_mat = diagm(ε)

    κ_current = κ
    κ_min = 0.001

    if verbose
        println("Iter     Free energy          ΔF          κ")
    end

    for iter in 1:maxiter
        F_old = F

        # gradients
        g_psi, Hmat = gradient_psi(psi, Hpsi, f, grid.dx)
        g_eta = gradient_eta(psi, Hmat, f, ε, kT, degeneracy)

        # preconditioned steepest descent directions
        Δ_psi = -apply_preconditioner(g_psi, precond_psi)
        Δ_eta = κ_current * (Hmat - diagm(ε))

        dot_psi = sum( dot(g_psi[:,n], Δ_psi[:,n]) for n=1:Nstates ) * grid.dx
        dot_eta = sum(g_eta .* Δ_eta)
        dot_current = dot_psi + dot_eta

        # reset CG if necessary
        if dot_current >= 0 || (iter > 1 && abs(dot_current/dot_old) < 1e-12)
            γ = 0.0
            X_psi .= 0.0
            X_eta .= 0.0
            dot_old = 0.0
            verbose && println("   CG reset at iter $iter")
        else
            γ = iter == 1 ? 0.0 : dot_current / dot_old
        end

        X_psi = Δ_psi + γ * X_psi
        X_eta = Δ_eta + γ * X_eta

        # ---------- line search ----------
        ξ_trial = 0.2
        psi_trial = psi + ξ_trial * X_psi
        eta_trial = eta_mat + ξ_trial * X_eta

        eig_trial = eigen(Symmetric(eta_trial))
        U_trial = eig_trial.vectors
        ε_trial = eig_trial.values

        psi_rot = psi_trial * U_trial
        orthonormalize!(psi_rot, grid.dx)

        μ_trial = find_mu(ε_trial, N_electrons, kT, degeneracy)
        f_trial = fermi.(ε_trial, μ_trial, kT, degeneracy)

        ρ_trial = calc_density(psi_rot, f_trial)
        V_H_trial = hartree_potential(ρ_trial, grid, a_soft_ee)
        V_tot_trial = V_ext + V_H_trial
        Hpsi_trial = apply_H(psi_rot, V_tot_trial, grid)
        Hmat_trial = dot_grid(psi_rot, Hpsi_trial, grid.dx)
        F_trial = calc_free_energy(psi_rot, f_trial, Hpsi_trial, Hmat_trial,
                               ρ_trial, V_H_trial, grid, V_ext, kT, degeneracy)

        dF0 = 2 * (dot_psi + dot_eta)
        A = 2 * (F_trial - F_old - dF0 * ξ_trial) / ξ_trial^2
        if A > 1e-8
            ξ_min = -dF0 / A
            ξ_min = clamp(ξ_min, -2ξ_trial, 2ξ_trial)
        else
            ξ_min = 0.1 * sign(dF0) * (dF0 > 0 ? -1 : 1)
        end

        # backtracking line search (Armijo condition) with μ monitoring
        α = 1.0
        for _ in 1:10
            psi_new = psi + α * ξ_min * X_psi
            eta_new = eta_mat + α * ξ_min * X_eta

            eig_new = eigen(Symmetric(eta_new))
            U_new = eig_new.vectors
            ε_new = eig_new.values

            psi_new_rot = psi_new * U_new
            orthonormalize!(psi_new_rot, grid.dx)

            μ_new = find_mu(ε_new, N_electrons, kT, degeneracy)
            f_new = fermi.(ε_new, μ_new, kT, degeneracy)

            ρ_new = calc_density(psi_new_rot, f_new)
            V_H_new = hartree_potential(ρ_new, grid, a_soft_ee)
            V_tot_new = V_ext + V_H_new
            Hpsi_new = apply_H(psi_new_rot, V_tot_new, grid)
            Hmat_new = dot_grid(psi_new_rot, Hpsi_new, grid.dx)
            F_new = calc_free_energy(psi_new_rot, f_new, Hpsi_new, Hmat_new,
                                 ρ_new, V_H_new, grid, V_ext, kT, degeneracy)

            # predicted change in μ (using w = f(1 - f/g))
            w = f .* (1.0 .- f / degeneracy)
            Sw = sum(w)
            if Sw > 0
                dμ_pred = sum( w .* diag(X_eta) ) / Sw * (α * ξ_min)
            else
                dμ_pred = 0.0
            end
            dμ_actual = μ_new - μ

            if abs(dμ_actual - dμ_pred) > 0.1 * abs(dμ_pred + 1e-12)
                α *= 0.5
                verbose && println("   μ deviation >10%, reducing step to α=$α")
            elseif F_new < F_old + 0.1 * α * ξ_min * dF0
                break
            else
                α *= 0.5
            end
        end

        # final update
        psi_new = psi + α * ξ_min * X_psi
        eta_new = eta_mat + α * ξ_min * X_eta

        eig_new = eigen(Symmetric(eta_new))
        U_new = eig_new.vectors
        ε_new = eig_new.values

        psi_new_rot = psi_new * U_new
        orthonormalize!(psi_new_rot, grid.dx)

        # rotate search directions
        X_psi = X_psi * U_new
        X_eta = U_new' * X_eta * U_new

        μ_new = find_mu(ε_new, N_electrons, kT, degeneracy)
        f_new = fermi.(ε_new, μ_new, kT, degeneracy)

        ρ_new = calc_density(psi_new_rot, f_new)
        V_H_new = hartree_potential(ρ_new, grid, a_soft_ee)
        V_tot_new = V_ext + V_H_new
        Hpsi_new = apply_H(psi_new_rot, V_tot_new, grid)
        Hmat_new = dot_grid(psi_new_rot, Hpsi_new, grid.dx)
        F_new = calc_free_energy(psi_new_rot, f_new, Hpsi_new, Hmat_new,
                             ρ_new, V_H_new, grid, V_ext, kT, degeneracy)

        if verbose
            println(iter, "    ", F_new, "    ", F_new - F_old, "    ", κ_current)
        end

        if abs(F_new - F_old) < tol
            verbose && println("Converged after $iter iterations")
            break
        end

        # adjust κ based on μ deviation
        w = f .* (1.0 .- f / degeneracy)
        Sw = sum(w)
        if Sw > 0
            dμ_pred = sum( w .* diag(X_eta) ) / Sw * (α * ξ_min)
            dμ_actual = μ_new - μ
            if abs(dμ_actual - dμ_pred) > 0.1 * abs(dμ_pred + 1e-12)
                κ_current = max(κ_current * 0.8, κ_min)
            else
                κ_current = min(κ_current * 1.05, 0.5)
            end
        end

        psi = psi_new_rot
        eta_mat = diagm(ε_new)
        ε = ε_new
        f = f_new
        μ = μ_new
        ρ = ρ_new
        V_H = V_H_new
        V_tot = V_tot_new
        Hpsi = Hpsi_new
        Hmat = Hmat_new
        F = F_new
        dot_old = dot_current
    end

    return psi, ε, f, μ, F
end

# ============================================================================
# Example runs
# ============================================================================
#=
function test_main(; psi=nothing, ebands=nothing)

    println("\n=== Spin‑unpolarised (g=2), two electrons ===")
    psi, ε, f, μ, F = minimize_dft(
        L=20.0, Ngrid=100,
        atoms=[-1.0, 1.0], Z=[1.0, 1.0],
        a_soft_ee=1.0, a_soft_en=1.0,
        N_electrons=2.0,
        kT=0.05,
        Nstates=8,
        degeneracy=2,
        κ=0.1,
        maxiter=500,
        tol=1e-6,
        verbose=true,
        seed=nothing,
        psi=psi,
        ε=ebands
    )
    serialize("psi.jldat", psi)
    serialize("ebands.jldat", ε)
    #
    println("\nFinal eigenvalues: ", ε)
    println("Occupations: ", f)
    println("Chemical potential: ", μ)
    println("Free energy: ", F)
end
=#