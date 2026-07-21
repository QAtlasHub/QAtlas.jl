# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: TFIM post-quench entanglement entropy S(ℓ, t).
#
# Implements the acceptance criteria of issue #144:
#
#  1.  t = 0 quench S(ℓ, 0) coincides with the equilibrium
#      VonNeumannEntropy of the *initial* Hamiltonian (back-compat
#      sanity).
#  2.  For h_0 = h_f (no quench) S(ℓ, t) is independent of t (basis
#      preserved by the evolution operator that diagonalises the GS).
#  3.  Quench h_0 = 2 → h_f = 1 (to critical), N = 32, ℓ = 8: S(ℓ, t)
#      grows monotonically over t ∈ [0, 4] and then saturates near a
#      volume-law value `≈ (c/6) ℓ + s_0`.
#  4.  A pinned-value sanity check at one (h_0, h_f, N, ℓ, t) point.
#  5.  Existing equilibrium API (`VonNeumannEntropy()`) still works
#      under the parametric type rewrite — back-compat for the
#      no-argument constructor.
#
# References:
#   * P. Calabrese, J. Cardy, J. Stat. Mech. P04010 (2005).
#   * I. Peschel, J. Phys. A 36, L205 (2003).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM QuenchEntanglementEntropy — back-compat sanity (t = 0)" begin
    # At t = 0 the post-quench reduced density matrix equals that of the
    # initial-Hamiltonian ground state; thus the quench S(ℓ, 0) must
    # match the *equilibrium* S(ℓ) of `initial`.
    N = 32
    for (h0, hf, ℓ) in ((2.0, 1.0, 8), (0.5, 1.0, 6), (1.0, 0.5, 10), (3.0, 0.7, 4))
        m_0 = TFIM(; J=1.0, h=h0)
        m_f = TFIM(; J=1.0, h=hf)
        S_eq_init = QAtlas.fetch(m_0, VonNeumannEntropy(), OBC(N); ℓ=ℓ)
        S_quench_t0 = QAtlas.fetch(
            m_f, QuenchEntanglementEntropy(), OBC(N); initial=m_0, ℓ=ℓ, t=0.0
        )
        @test S_quench_t0 ≈ S_eq_init atol = 1e-12
    end
end

@testset "TFIM QuenchEntanglementEntropy — h_0 = h_f leaves S(ℓ, t) invariant" begin
    # The ground state is an eigenstate of H, so |Ψ(t)⟩ = exp(-i E_0 t) |Ψ_0⟩
    # and the reduced density matrix is unchanged.  Numerically, the
    # Majorana congruence Σ → R Σ R^T must reduce to a similarity that
    # preserves the spectrum of i Σ_A.
    N = 32
    ℓ = 8
    for h in (0.5, 1.0, 2.0)
        m = TFIM(; J=1.0, h=h)
        S_t0 = QAtlas.fetch(m, QuenchEntanglementEntropy(), OBC(N); initial=m, ℓ=ℓ, t=0.0)
        for t in (0.7, 1.5, 3.0, 5.0)
            S_t = QAtlas.fetch(m, QuenchEntanglementEntropy(), OBC(N); initial=m, ℓ=ℓ, t=t)
            @test S_t ≈ S_t0 atol = 1e-9
        end
    end
end

@testset "TFIM QuenchEntanglementEntropy — monotone linear growth then saturation" begin
    # Calabrese–Cardy quasi-particle picture: for a quench from an
    # initial gapped TFIM into the gapless TFIM (h_f = J = 1), S(ℓ, t)
    # grows ≈ linearly until t* = ℓ / (2 v_E) and saturates at a
    # volume-law value.  We test the qualitative shape:
    #
    #   * monotone non-decreasing on t ∈ [0, 4] (the linear-growth
    #     window),
    #   * saturation: S(t = 8) and S(t = 12) are both substantially
    #     above the t = 0 value and within ~10 % of each other (finite-N
    #     revival oscillation),
    #   * the saturation magnitude is order `(c/6) ℓ ≈ 0.66` (with
    #     `c = 1/2`, `ℓ = 8`) plus the non-universal initial-state
    #     constant; it falls in [1.0, 3.0].

    m_0 = TFIM(; J=1.0, h=2.0)
    m_f = TFIM(; J=1.0, h=1.0)   # critical, c = 1/2
    N = 32
    ℓ = 8

    ts_grow = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0)
    Ss_grow = [
        QAtlas.fetch(m_f, QuenchEntanglementEntropy(), OBC(N); initial=m_0, ℓ=ℓ, t=t) for
        t in ts_grow
    ]
    # Strict monotonic increase over the linear-growth window.
    for k in 2:length(Ss_grow)
        @test Ss_grow[k] > Ss_grow[k - 1]
    end

    S_sat_a = QAtlas.fetch(
        m_f, QuenchEntanglementEntropy(), OBC(N); initial=m_0, ℓ=ℓ, t=8.0
    )
    S_sat_b = QAtlas.fetch(
        m_f, QuenchEntanglementEntropy(), OBC(N); initial=m_0, ℓ=ℓ, t=12.0
    )
    # Saturation: well above the small initial value and within finite-N
    # revival window of each other.
    @test S_sat_a > Ss_grow[end]                 # still above S(t = 4)
    @test S_sat_b > Ss_grow[1] + 1.0             # >> S(t = 0)
    @test abs(S_sat_a - S_sat_b) / S_sat_a < 0.1  # finite-N revivals < 10 %
    # Order-of-magnitude check on the volume-law plateau.
    @test 1.0 < S_sat_a < 3.0

    # Coarse linear-growth slope estimate (least squares) on (0, 0.5, …, 3.0):
    # (c/3) v_E ≈ (1/6) · 2 ≈ 0.33 with v_E = 2J = 2 at h_f = J.  We
    # only assert that the slope sits in a generous physically motivated
    # window; the precise value depends on the (non-universal) initial
    # state via the Bogoliubov occupations.
    n = 7
    xs = collect(ts_grow[1:n])
    ys = Ss_grow[1:n]
    x̄ = sum(xs) / n
    ȳ = sum(ys) / n
    slope = sum((xs[i] - x̄) * (ys[i] - ȳ) for i in 1:n) / sum((xs[i] - x̄)^2 for i in 1:n)
    @test 0.2 < slope < 0.6
end

@testset "TFIM QuenchEntanglementEntropy — pinned value (regression)" begin
    # Pin one numerical value computed at the time the implementation
    # was written.  The atol is chosen comfortably above the round-off
    # ceiling of the matrix exponential + Hermitian eigendecomposition
    # (~1e-12 relative).
    m_0 = TFIM(; J=1.0, h=2.0)
    m_f = TFIM(; J=1.0, h=1.0)
    S = QAtlas.fetch(m_f, QuenchEntanglementEntropy(), OBC(32); initial=m_0, ℓ=8, t=1.0)
    @test S ≈ 0.5338124210270956 atol = 1e-10
end

@testset "TFIM QuenchEntanglementEntropy — input validation" begin
    m_0 = TFIM(; J=1.0, h=2.0)
    m_f = TFIM(; J=1.0, h=1.0)
    @test_throws ArgumentError QAtlas.fetch(
        m_f, QuenchEntanglementEntropy(), OBC(8); initial=m_0, ℓ=0, t=1.0
    )
    @test_throws ArgumentError QAtlas.fetch(
        m_f, QuenchEntanglementEntropy(), OBC(8); initial=m_0, ℓ=8, t=1.0
    )
end

@testset "TFIM VonNeumannEntropy — equilibrium back-compat after the mode split" begin
    # #734: `VonNeumannEntropy` is now AbstractQAtlas's singleton and the
    # post-quench branch is the separate `QuenchEntanglementEntropy`.  Guard the
    # ergonomic `VonNeumannEntropy()` API used throughout the test suite: it must
    # still mean the equilibrium entropy and still deliver the Peschel value.
    N = 16
    m = TFIM(; J=1.0, h=1.0)
    @test VonNeumannEntropy() isa VonNeumannEntropy
    @test QuenchEntanglementEntropy() isa QuenchEntanglementEntropy
    # The two branches are now distinct types, not two instantiations of one.
    @test VonNeumannEntropy !== QuenchEntanglementEntropy
    @test QuenchEntanglementEntropy <: QAtlas.AbstractEntanglementMeasure
    # The mode-keyed constructor is gone with the parametric type (no shim: a
    # `VonNeumannEntropy(::Symbol)` method here would be piracy on AbstractQAtlas).
    @test_throws MethodError VonNeumannEntropy(:quench)

    S = QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); ℓ=8)
    @test S > 0  # non-trivial entanglement at criticality
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "TFIM quench entanglement — verification cards" begin
    # t = 0: the quench entanglement entropy equals the entanglement of
    # the initial ground state (independent route: equilibrium S_vN).
    let m0 = TFIM(; J=1.0, h=2.0), mf = TFIM(; J=1.0, h=0.5), N = 8, ℓ = 4
        verify(
            mf,
            QuenchEntanglementEntropy(),
            OBC(N);
            route=:limiting_case,
            fetch_kw=(; initial=m0, ℓ=ℓ, t=0.0),
            independent=QAtlas.fetch(m0, VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=Inf),
            agree_within=1e-8,
            refs=["t=0: quench EE equals the initial ground-state entanglement"],
        )
    end
end
