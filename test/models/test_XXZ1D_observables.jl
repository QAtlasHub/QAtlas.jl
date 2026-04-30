using QAtlas, Test, LinearAlgebra, ForwardDiff

# =============================================================================
# Tests for the XXZ1D dense-ED finite-T observables added in Tier-1 of the
# observable-coverage parity work (TFIM ↔ XXZ1D / Heisenberg1D).  Layered
# the same way as `test_TFIM_thermal.jl`:
#
#   1. T → 0:  ⟨A⟩_β at large β ≈ ground-state value.
#   2. T → ∞:  high-temperature paramagnet limits.
#   3. ED comparison at small N (independent kron / direct ρ build).
#   4. SU(2) symmetry checks at the Heisenberg point Δ = 1.
#   5. Entanglement: subsystem entropy bounds + ground-state SVD baseline.
# =============================================================================

# --- Independent dense-ED helpers (do NOT use QAtlas internals) -------------

const _SX_dir = ComplexF64[0 1; 1 0]
const _SY_dir = ComplexF64[0 -im; im 0]
const _SZ_dir = ComplexF64[1 0; 0 -1]
const _I2_dir = ComplexF64[1 0; 0 1]

function _direct_pauli_at(N::Int, i::Int, σ::Matrix{ComplexF64})
    ops = [j == i ? σ : _I2_dir for j in 1:N]
    return reduce(kron, ops)
end

function _direct_pair_at(N::Int, i::Int, j::Int, σ::Matrix{ComplexF64})
    ops = [(k == i || k == j) ? σ : _I2_dir for k in 1:N]
    return reduce(kron, ops)
end

function _direct_xxz_H(N::Int, J::Float64, Δ::Float64)
    D = 2^N
    H = zeros(ComplexF64, D, D)
    for i in 1:(N - 1)
        H .+= (J / 4) .* _direct_pair_at(N, i, i + 1, _SX_dir)
        H .+= (J / 4) .* _direct_pair_at(N, i, i + 1, _SY_dir)
        H .+= (J * Δ / 4) .* _direct_pair_at(N, i, i + 1, _SZ_dir)
    end
    return H
end

function _direct_thermal_state(H::AbstractMatrix, β::Real)
    F = eigen(Hermitian(H))
    emin = minimum(F.values)
    ws = exp.(-β .* (F.values .- emin))
    Z = sum(ws)
    ρ = F.vectors * Diagonal(ComplexF64.(ws ./ Z)) * F.vectors'
    return ρ, F
end

# ───────────────────────────────────────────────────────────────────────
# Layer 2: high-temperature paramagnet limits
# ───────────────────────────────────────────────────────────────────────
@testset "XXZ1D T → ∞ paramagnet limits" begin
    β_high = 1e-4
    for Δ in (-0.5, 0.0, 0.7, 1.0), N in (4, 5)
        m = XXZ1D(; J=1.0, Δ=Δ)

        # All bulk magnetisations vanish at infinite T (Tr σᵅ = 0).
        @test abs(QAtlas.fetch(m, MagnetizationX(), OBC(N); beta=β_high)) < 5e-5
        @test abs(QAtlas.fetch(m, MagnetizationY(), OBC(N); beta=β_high)) < 1e-12
        @test abs(QAtlas.fetch(m, MagnetizationZ(), OBC(N); beta=β_high)) < 1e-12

        # Each-site σᵅ has variance ⟨σᵅ²⟩ = 1, off-site connected
        # contributions vanish at β → 0, so χ_αα → β → 0.
        @test QAtlas.fetch(m, SusceptibilityXX(), OBC(N); beta=β_high) ≈ β_high atol = 5e-6
        @test QAtlas.fetch(m, SusceptibilityYY(), OBC(N); beta=β_high) ≈ β_high atol = 5e-6
        @test QAtlas.fetch(m, SusceptibilityZZ(), OBC(N); beta=β_high) ≈ β_high atol = 5e-6
    end
end

# ───────────────────────────────────────────────────────────────────────
# Layer 1: T → 0 — large-β thermal averages match ground state
# ───────────────────────────────────────────────────────────────────────
@testset "XXZ1D T → 0 collapses to ground state observables" begin
    β = 80.0
    for (J, Δ, N) in ((1.0, 0.5, 5), (1.3, 1.0, 5))
        m = XXZ1D(; J=J, Δ=Δ)

        # MagnetizationY at β = ∞ is exactly zero from any Hermitian-real H.
        my = QAtlas.fetch(m, MagnetizationY(), OBC(N); beta=β)
        @test abs(my) < 1e-10

        # Mz at GS: U(1) symmetry ⇒ ground state has definite Sᶻ; even N
        # has Sᶻ = 0 sector ground state, odd N has Sᶻ = ±1/2 with twofold
        # ground-state degeneracy whose symmetric average is again zero.
        # The Boltzmann sum picks the degenerate manifold so Mz averages to 0.
        mz = QAtlas.fetch(m, MagnetizationZ(), OBC(N); beta=β)
        @test abs(mz) < 1e-8
    end
end

# ───────────────────────────────────────────────────────────────────────
# Layer 3: independent ED cross-check at small N
# ───────────────────────────────────────────────────────────────────────
@testset "XXZ1D ED cross-check (N=3, 4)" begin
    for (J, Δ, N, β) in ((1.0, 0.5, 3, 1.0), (1.3, 0.7, 4, 0.7), (0.8, 1.0, 4, 2.0))
        m = XXZ1D(; J=J, Δ=Δ)
        H = _direct_xxz_H(N, J, Δ)
        ρ, _ = _direct_thermal_state(H, β)

        Mx = sum(_direct_pauli_at(N, i, _SX_dir) for i in 1:N)
        My = sum(_direct_pauli_at(N, i, _SY_dir) for i in 1:N)
        Mz = sum(_direct_pauli_at(N, i, _SZ_dir) for i in 1:N)
        mx_ed = real(tr(ρ * Mx)) / N
        my_ed = real(tr(ρ * My)) / N
        mz_ed = real(tr(ρ * Mz)) / N

        @test QAtlas.fetch(m, MagnetizationX(), OBC(N); beta=β) ≈ mx_ed atol = 1e-10
        @test QAtlas.fetch(m, MagnetizationY(), OBC(N); beta=β) ≈ my_ed atol = 1e-10
        @test QAtlas.fetch(m, MagnetizationZ(), OBC(N); beta=β) ≈ mz_ed atol = 1e-10

        # Susceptibilities χ_αα = β · Var(M_α) / N
        χxx_ed = β * (real(tr(ρ * (Mx * Mx))) - real(tr(ρ * Mx))^2) / N
        χyy_ed = β * (real(tr(ρ * (My * My))) - real(tr(ρ * My))^2) / N
        χzz_ed = β * (real(tr(ρ * (Mz * Mz))) - real(tr(ρ * Mz))^2) / N
        @test QAtlas.fetch(m, SusceptibilityXX(), OBC(N); beta=β) ≈ χxx_ed atol = 1e-10
        @test QAtlas.fetch(m, SusceptibilityYY(), OBC(N); beta=β) ≈ χyy_ed atol = 1e-10
        @test QAtlas.fetch(m, SusceptibilityZZ(), OBC(N); beta=β) ≈ χzz_ed atol = 1e-10

        # Two-point correlators (all axes, both modes)
        for (qty, σ) in
            ((XXCorrelation, _SX_dir), (YYCorrelation, _SY_dir), (ZZCorrelation, _SZ_dir))
            for i in 1:N, j in i:N
                if i == j
                    continue  # σᵅ² = I, skipped
                end
                Op = _direct_pair_at(N, i, j, σ)
                Oi = _direct_pauli_at(N, i, σ)
                Oj = _direct_pauli_at(N, j, σ)
                c2 = real(tr(ρ * Op))
                ci = real(tr(ρ * Oi))
                cj = real(tr(ρ * Oj))
                v_static = QAtlas.fetch(m, qty(; mode=:static), OBC(N); beta=β, i=i, j=j)
                v_conn = QAtlas.fetch(m, qty(; mode=:connected), OBC(N); beta=β, i=i, j=j)
                @test v_static ≈ c2 atol = 1e-10
                @test v_conn ≈ c2 - ci * cj atol = 1e-10
            end
        end

        # EnergyLocal: bonds split symmetrically must sum to ⟨H⟩
        ε = QAtlas.fetch(m, EnergyLocal(), OBC(N); beta=β)
        E_total = QAtlas.fetch(m, Energy(:total), OBC(N); beta=β)
        @test sum(ε) ≈ E_total atol = 1e-10

        # Local magnetisations (length N, near-zero by U(1) symmetry)
        mxv = QAtlas.fetch(m, MagnetizationXLocal(), OBC(N); beta=β)
        myv = QAtlas.fetch(m, MagnetizationYLocal(), OBC(N); beta=β)
        mzv = QAtlas.fetch(m, MagnetizationZLocal(), OBC(N); beta=β)
        @test length(mxv) == N
        @test length(myv) == N
        @test length(mzv) == N
        @test all(abs.(mxv) .< 1e-10)
        @test all(abs.(myv) .< 1e-10)
        @test all(abs.(mzv) .< 1e-8)
    end
end

# ───────────────────────────────────────────────────────────────────────
# SU(2) symmetry: at Δ = 1 (Heisenberg) χ_xx = χ_yy = χ_zz
# ───────────────────────────────────────────────────────────────────────
@testset "XXZ1D Heisenberg point: SU(2)-symmetric susceptibilities" begin
    m = XXZ1D(; J=1.0, Δ=1.0)
    for β in (0.5, 1.5, 3.0), N in (4, 5)
        χxx = QAtlas.fetch(m, SusceptibilityXX(), OBC(N); beta=β)
        χyy = QAtlas.fetch(m, SusceptibilityYY(), OBC(N); beta=β)
        χzz = QAtlas.fetch(m, SusceptibilityZZ(), OBC(N); beta=β)
        @test χxx ≈ χyy atol = 1e-10
        @test χyy ≈ χzz atol = 1e-10
    end
end

# ───────────────────────────────────────────────────────────────────────
# Mass gap
# ───────────────────────────────────────────────────────────────────────
@testset "XXZ1D MassGap" begin
    # OBC: matches direct spectrum diff
    for (J, Δ, N) in ((1.0, 0.5, 5), (1.3, 1.0, 4), (0.8, -0.3, 4))
        m = XXZ1D(; J=J, Δ=Δ)
        H = QAtlas._xxz1d_hamiltonian_matrix(m, N)
        evals = sort(real.(eigvals(Hermitian(H))))
        gap_direct = evals[2] - evals[1]
        @test QAtlas.fetch(m, MassGap(), OBC(N)) ≈ gap_direct atol = 1e-12
    end

    # Infinite: critical regime → 0; gapped regime → NaN with warning
    @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.0), MassGap(), Infinite()) == 0.0
    @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.5), MassGap(), Infinite()) == 0.0
    @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), MassGap(), Infinite()) == 0.0
    val = @test_logs (:warn, r"gapped-regime") QAtlas.fetch(
        XXZ1D(; J=1.0, Δ=1.5), MassGap(), Infinite()
    )
    @test isnan(val)
end

# ───────────────────────────────────────────────────────────────────────
# Entanglement entropies: bounds + GS SVD reference
# ───────────────────────────────────────────────────────────────────────
@testset "XXZ1D entanglement entropies — bounds and SVD reference" begin
    # Ground-state EE matches the direct SVD of the GS vector |ψ⟩.
    for (J, Δ, N, ℓ) in ((1.0, 0.5, 4, 2), (1.0, 1.0, 5, 2), (1.0, 0.0, 5, 2))
        m = XXZ1D(; J=J, Δ=Δ)
        H = QAtlas._xxz1d_hamiltonian_matrix(m, N)
        F = eigen(Hermitian(H))
        ψ = F.vectors[:, 1]
        Ψ = reshape(ψ, (2^ℓ, 2^(N - ℓ)))
        svals = svdvals(Ψ)
        S_direct = -sum(s -> begin
            p = s^2
            p > 1e-15 ? p * log(p) : 0.0
        end, svals)

        S_qa = QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=Inf)
        @test S_qa ≈ S_direct atol = 1e-9
        # 0 ≤ S ≤ ℓ log 2
        @test 0.0 ≤ S_qa ≤ ℓ * log(2) + 1e-10

        # Renyi(2) ≤ vN ≤ ℓ log 2.
        S2 = QAtlas.fetch(m, RenyiEntropy(2), OBC(N); ℓ=ℓ, beta=Inf)
        @test 0.0 ≤ S2 ≤ S_qa + 1e-10

        # Renyi(α → 0): S_0 = log(rank ρ_A) ≤ ℓ log 2.
        S_half = QAtlas.fetch(m, RenyiEntropy(0.5), OBC(N); ℓ=ℓ, beta=Inf)
        @test S_half ≥ S_qa - 1e-10
    end

    # Thermal mixed state: S(β=Inf) should match S(β=80) up to gap-suppressed
    # corrections, and S(β=0+) → ℓ log 2 (maximally mixed).
    let m = XXZ1D(; J=1.0, Δ=0.5), N = 4, ℓ = 2
        S_inf = QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=Inf)
        S_lo = QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=80.0)
        @test S_inf ≈ S_lo atol = 1e-6

        S_hi = QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=1e-4)
        @test S_hi ≈ ℓ * log(2) atol = 1e-3
    end
end

# ───────────────────────────────────────────────────────────────────────
# AutoDiff sanity: c_v = -β² ∂ε/∂β  on XXZ1D OBC
# ───────────────────────────────────────────────────────────────────────
@testset "XXZ1D specific heat from AutoDiff" begin
    for (Δ, N) in ((-0.5, 4), (0.0, 5), (0.7, 4), (1.0, 5)), β in (0.5, 1.5)
        m = XXZ1D(; J=1.0, Δ=Δ)
        ε = b -> QAtlas.fetch(m, Energy(:per_site), OBC(N); beta=b)
        c_ad = -β^2 * ForwardDiff.derivative(ε, β)
        c_var = QAtlas.fetch(m, SpecificHeat(), OBC(N); beta=β)
        @test c_ad ≈ c_var atol = 1e-9 rtol = 1e-9
    end
end
