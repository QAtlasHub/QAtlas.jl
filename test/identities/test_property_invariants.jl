using Test
using QAtlas
using Random: Random, seed!

# Property-based invariant tests — randomly sample (model parameters,
# β, N) from a deterministic seed and assert universal invariants of
# the corresponding fetch outputs.  The harness's strength is *coverage*:
# each randomly-drawn slice spot-checks a different (J, h, Δ, β, N)
# combination, accumulating confidence that the implementations satisfy
# basic physical bounds across the phase diagram, not just at the
# canonical points hard-coded in the per-model identity files.
#
# The seed is fixed (`Random.seed!(0xQAtlas)`) so the sweep is
# reproducible across CI runs and reproduces the same failure on
# regression.  Increase the iteration count if a more aggressive
# fuzz is needed; default `n_samples = 20` keeps wall time below
# ~5 s per model.

const _PROP_SEED = 0xCA7_A7_A5  # arbitrary but fixed
const _N_SAMPLES = 20

# Helpers
_inrange(x, a, b) = (a ≤ x ≤ b) || isapprox(x, a; atol=1e-10) || isapprox(x, b; atol=1e-10)

# ────────────────────────────────────────────────────────────────────
# (1) Thermodynamic stability:  c_v ≥ 0 and s ≥ 0  (any model, any β)
# ────────────────────────────────────────────────────────────────────

@testset "TFIM invariants — random (h, β, N) sweep" begin
    seed!(_PROP_SEED)
    for k in 1:_N_SAMPLES
        h = 0.2 + 2.6 * rand()       # h ∈ [0.2, 2.8]
        J = 1.0
        β = exp(-1 + 3 * rand())     # β ∈ [exp(-1), exp(2)] ≈ [0.37, 7.4]
        N = rand((4, 6, 8))
        model = TFIM(; J=J, h=h)
        bc = OBC(N)

        ε = QAtlas.fetch(model, Energy(:per_site), bc; beta=β)
        f = QAtlas.fetch(model, FreeEnergy(), bc; beta=β)
        s = QAtlas.fetch(model, ThermalEntropy(), bc; beta=β)
        c = QAtlas.fetch(model, SpecificHeat(), bc; beta=β)
        m_x = QAtlas.fetch(model, MagnetizationX(), bc; beta=β)
        χ_xx = QAtlas.fetch(model, SusceptibilityXX(), bc; beta=β)

        # Thermodynamic stability
        @test c ≥ -1e-12
        @test s ≥ -1e-12

        # Operator bounds: |⟨σˣ⟩| ≤ 1
        @test abs(m_x) ≤ 1 + 1e-12

        # Variance non-negativity (χ_xx OBC is β·Var(M_x)/N ≥ 0).
        @test χ_xx ≥ -1e-12

        # Energy bounds: ε ≥ E_GS / N (OBC) and ε ≤ 0 at β > 0 since
        # H = -J σᶻσᶻ - h σˣ has Tr H = 0 → ε(β=0) = 0 → ε(β>0) ≤ 0.
        @test ε ≤ 1e-12

        # Free energy ≤ Energy: f = ε - T·s ≤ ε since s ≥ 0.
        @test f ≤ ε + 1e-12
    end
end

# ────────────────────────────────────────────────────────────────────
# (2) Cooling monotonicity:  ε(β) decreases in β
# ────────────────────────────────────────────────────────────────────

@testset "TFIM cooling: ε(β) monotonically decreases in β" begin
    seed!(_PROP_SEED + 1)
    for k in 1:10
        h = 0.3 + 2.0 * rand()
        N = rand((6, 8))
        model = TFIM(; J=1.0, h=h)
        bc = OBC(N)
        βs = sort(rand(4) .* 5 .+ 0.1)  # 4 ascending β's
        εs = [QAtlas.fetch(model, Energy(:per_site), bc; beta=β) for β in βs]
        @test issorted(εs; rev=true)
    end
end

# ────────────────────────────────────────────────────────────────────
# (3) High-T / low-T entropy limits
# ────────────────────────────────────────────────────────────────────

@testset "TFIM entropy limits: s(β → 0) → log 2,  s(β → ∞) → 0" begin
    seed!(_PROP_SEED + 2)
    for k in 1:8
        h = 0.3 + 2.0 * rand()
        N = rand((6, 8))
        model = TFIM(; J=1.0, h=h)
        bc = OBC(N)
        # High-T limit: at β → 0 the thermal state is `I/2^N`, so per-
        # site entropy → log 2.
        s_hi = QAtlas.fetch(model, ThermalEntropy(), bc; beta=1e-4)
        @test s_hi ≈ log(2) atol = 1e-3
        # Low-T limit: at β → ∞ the GS dominates → per-site entropy → 0.
        # OBC TFIM at finite N has an edge-mode degeneracy that gives
        # ~log(2)/N residual entropy near criticality (`|h − J|` small),
        # so we only assert the low-T value is well below the high-T
        # log(2) ceiling, not literally zero.
        s_lo = QAtlas.fetch(model, ThermalEntropy(), bc; beta=80.0)
        @test 0 ≤ s_lo < 0.2
    end
end

# ────────────────────────────────────────────────────────────────────
# (4) XXZ1D invariants  (random Δ ∈ (-1, 1), β, N)
# ────────────────────────────────────────────────────────────────────

@testset "XXZ1D invariants — random (Δ, β, N) sweep" begin
    seed!(_PROP_SEED + 10)
    for k in 1:_N_SAMPLES
        Δ = -0.9 + 1.8 * rand()      # critical regime
        β = exp(-1 + 2 * rand())
        N = rand((4, 6, 8))
        model = XXZ1D(; J=1.0, Δ=Δ)
        bc = OBC(N)

        f = QAtlas.fetch(model, FreeEnergy(), bc; beta=β)
        s = QAtlas.fetch(model, ThermalEntropy(), bc; beta=β)
        c = QAtlas.fetch(model, SpecificHeat(), bc; beta=β)
        m_x = QAtlas.fetch(model, MagnetizationX(), bc; beta=β)
        m_y = QAtlas.fetch(model, MagnetizationY(), bc; beta=β)
        m_z = QAtlas.fetch(model, MagnetizationZ(), bc; beta=β)
        χ_xx = QAtlas.fetch(model, SusceptibilityXX(), bc; beta=β)
        χ_yy = QAtlas.fetch(model, SusceptibilityYY(), bc; beta=β)
        χ_zz = QAtlas.fetch(model, SusceptibilityZZ(), bc; beta=β)

        @test s ≥ -1e-12
        @test c ≥ -1e-12
        @test χ_xx ≥ -1e-12
        @test χ_yy ≥ -1e-12
        @test χ_zz ≥ -1e-12
        # |⟨σᵅ⟩| ≤ 1 per site
        @test abs(m_x) ≤ 1 + 1e-12
        @test abs(m_y) ≤ 1 + 1e-12
        @test abs(m_z) ≤ 1 + 1e-12
        # XXZ1D in canonical ensemble has m_α = 0 by construction
        # (no symmetry-breaking field, U(1) × Z₂ → m_x, m_y, m_z all
        # vanish at finite N at any β).
        @test abs(m_x) < 1e-10
        @test abs(m_y) < 1e-10
        @test abs(m_z) < 1e-10
    end
end

# ────────────────────────────────────────────────────────────────────
# (5) Entanglement entropy bounds:  0 ≤ S_vN ≤ ℓ log 2
# ────────────────────────────────────────────────────────────────────

@testset "TFIM entanglement bounds: 0 ≤ S_vN ≤ ℓ log 2" begin
    seed!(_PROP_SEED + 20)
    for k in 1:10
        h = 0.3 + 2.0 * rand()
        N = rand((8, 12))
        ℓ = rand(2:(N - 2))
        β = exp(-1 + 2 * rand())
        S = QAtlas.fetch(TFIM(; J=1.0, h=h), VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=β)
        @test S ≥ -1e-12
        # ℓ log 2 is the maximally-mixed bound on the half-chain reduced
        # density matrix.
        @test S ≤ ℓ * log(2) + 1e-10
    end
end

@testset "TFIM Renyi bounds: 0 ≤ S_α ≤ ℓ log 2,  S_α ≤ S_vN for α > 1" begin
    seed!(_PROP_SEED + 21)
    for k in 1:10
        h = 0.3 + 2.0 * rand()
        N = rand((8, 12))
        ℓ = rand(2:(N - 2))
        # Random α in (1, 5) — all Renyi entropies for α > 1 are
        # bounded above by S_vN.
        α = 1.1 + 4.0 * rand()
        S_vn = QAtlas.fetch(TFIM(; J=1.0, h=h), VonNeumannEntropy(), OBC(N); ℓ=ℓ)
        S_α = QAtlas.fetch(TFIM(; J=1.0, h=h), RenyiEntropy(α), OBC(N); ℓ=ℓ)
        @test S_α ≥ -1e-12
        @test S_α ≤ ℓ * log(2) + 1e-10
        @test S_α ≤ S_vn + 1e-10
    end
end

# ────────────────────────────────────────────────────────────────────
# (6) S1Heisenberg1D invariants  (spin-1, dim = 3)
# ────────────────────────────────────────────────────────────────────

@testset "S1Heisenberg1D invariants — random (J, β, N) sweep" begin
    seed!(_PROP_SEED + 30)
    # 3^N Hilbert space; cap N at 4 to keep the sweep cheap.
    for k in 1:10
        J = 0.3 + 1.2 * rand()
        β = exp(-1 + 2 * rand())
        N = rand((3, 4))
        model = S1Heisenberg1D(; J=J)
        bc = OBC(N)

        s = QAtlas.fetch(model, ThermalEntropy(), bc; beta=β)
        c = QAtlas.fetch(model, SpecificHeat(), bc; beta=β)

        @test s ≥ -1e-12
        @test c ≥ -1e-12
        # spin-1: |⟨S^α⟩| ≤ 1 (eigenvalues are 0, ±1)
        m_x = QAtlas.fetch(model, MagnetizationX(), bc; beta=β)
        m_z = QAtlas.fetch(model, MagnetizationZ(), bc; beta=β)
        @test abs(m_x) ≤ 1 + 1e-12
        @test abs(m_z) ≤ 1 + 1e-12

        # SU(2) symmetry: m_α = 0 in canonical ensemble at any (J, β).
        @test abs(m_x) < 1e-10
        @test abs(m_z) < 1e-10
    end
end

@testset "S1Heisenberg1D high-T entropy: s(β → 0) → log 3 (spin-1 dim)" begin
    seed!(_PROP_SEED + 31)
    for k in 1:5
        J = 0.5 + rand()
        N = rand((3, 4))
        s = QAtlas.fetch(S1Heisenberg1D(; J=J), ThermalEntropy(), OBC(N); beta=1e-4)
        @test s ≈ log(3) atol = 1e-3
    end
end

# ────────────────────────────────────────────────────────────────────
# (7) Static correlator Cauchy-Schwarz: |⟨σᵅᵢ σᵅⱼ⟩| ≤ 1
# ────────────────────────────────────────────────────────────────────

@testset "Static σᵅσᵅ correlators: |⟨σᵅ_i σᵅ_j⟩| ≤ 1" begin
    seed!(_PROP_SEED + 40)
    for k in 1:_N_SAMPLES
        h = 0.3 + 2.0 * rand()
        β = exp(-1 + 2 * rand())
        N = rand((6, 8))
        i = rand(2:(N - 1))
        j = rand(2:(N - 1))
        c_zz = QAtlas.fetch(
            TFIM(; J=1.0, h=h), ZZCorrelation{:static}(), OBC(N); beta=β, i=i, j=j
        )
        c_xx = QAtlas.fetch(
            TFIM(; J=1.0, h=h), XXCorrelation{:static}(), OBC(N); beta=β, i=i, j=j
        )
        @test abs(c_zz) ≤ 1 + 1e-10
        @test abs(c_xx) ≤ 1 + 1e-10
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "property invariants — verification cards" begin
    # Gibbs identity s = β(ε − f): the entropy is cross-checked against
    # the independently-implemented Energy and FreeEnergy via the
    # thermodynamic sum rule (three separate code paths must agree).
    let m = TFIM(; J=1.0, h=0.5), β = 1.3, N = 6
        ε = QAtlas.fetch(m, Energy(:per_site), OBC(N); beta=β)
        f = QAtlas.fetch(m, FreeEnergy(), OBC(N); beta=β)
        verify(
            m,
            ThermalEntropy(),
            OBC(N);
            route=:sum_rule,
            fetch_kw=(; beta=β),
            independent=β * (ε - f),
            agree_within=1e-7,
            refs=["Gibbs identity s = β(ε − f) across Energy/FreeEnergy/Entropy paths"],
        )
    end
end
