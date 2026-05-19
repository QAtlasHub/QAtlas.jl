# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: TFIM Loschmidt echo + DQPT rate function (issue #143).
#
# Validates `LoschmidtEcho{:amplitude}` (OBC), `LoschmidtEcho{:rate}` (OBC and
# Infinite) for sudden quenches `H_0 = TFIM(J, h_0) → H_f = TFIM(J, h_f)`.
#
# Reference: Heyl–Polkovnikov–Kehrein, PRL 110, 135704 (2013); Heyl, Rep. Prog.
# Phys. 81, 054001 (2018).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

# ─── small helper: locate k* and t_c for an Infinite paramag→ferromag quench ───

# Bogoliubov 2θ_k(h) = atan2(J sin k, h - J cos k) — same convention as the
# implementation in src/models/quantum/TFIM/TFIM_loschmidt.jl.
_two_theta(J, h, k) = atan(J * sin(k), h - J * cos(k))
_cos_two_dtheta(J, h0, hf, k) = cos(_two_theta(J, h0, k) - _two_theta(J, hf, k))
_Λ(J, h, k) = 2 * sqrt(J^2 + h^2 - 2 * J * h * cos(k))

# Bisect for k* ∈ (a, b) where cos(2 Δθ_{k*}) = 0; assumes a sign change in
# the bracket.
function _find_kstar(J, h0, hf; a=1e-6, b=π - 1e-6, tol=1e-12, maxit=200)
    fa = _cos_two_dtheta(J, h0, hf, a)
    fb = _cos_two_dtheta(J, h0, hf, b)
    @assert sign(fa) != sign(fb) "no sign change of cos(2Δθ) on (a, b); h_0,h_f must straddle the QCP"
    for _ in 1:maxit
        m = 0.5 * (a + b)
        fm = _cos_two_dtheta(J, h0, hf, m)
        if abs(fm) < tol || (b - a) < tol
            return m
        end
        if sign(fm) == sign(fa)
            a, fa = m, fm
        else
            b, fb = m, fm
        end
    end
    return 0.5 * (a + b)
end

@testset "TFIM Loschmidt echo (issue #143)" begin
    J = 1.0

    # ─── L(0) = 1, λ(0) = 0 ────────────────────────────────────────────────
    @testset "t = 0 invariants" begin
        m_0 = TFIM(; J=J, h=2.0)
        m_f = TFIM(; J=J, h=0.5)
        N = 32
        @test QAtlas.fetch(m_f, LoschmidtEcho(:amplitude), OBC(N); initial=m_0, t=0.0) ≈ 1.0 atol =
            1e-12
        @test QAtlas.fetch(m_f, LoschmidtRateFunction(), OBC(N); initial=m_0, t=0.0) ≈ 0.0 atol =
            1e-12
        @test QAtlas.fetch(m_f, LoschmidtRateFunction(), Infinite(); initial=m_0, t=0.0) ≈
            0.0 atol = 1e-12
    end

    # ─── No-quench identity: h_0 = h_f ⇒ L(t) = 1, λ(t) = 0 for all t ─────
    @testset "no-quench identity (h_0 = h_f)" begin
        m = TFIM(; J=J, h=1.3)
        N = 24
        for t in (0.1, 0.7, 2.5, 5.0)
            @test QAtlas.fetch(m, LoschmidtEcho(:amplitude), OBC(N); initial=m, t=t) ≈ 1.0 atol =
                1e-10
            @test QAtlas.fetch(m, LoschmidtRateFunction(), OBC(N); initial=m, t=t) ≈ 0.0 atol =
                1e-10
            @test QAtlas.fetch(m, LoschmidtRateFunction(), Infinite(); initial=m, t=t) ≈ 0.0 atol =
                1e-10
        end
    end

    # ─── L(t) ∈ [0, 1], λ(t) ≥ 0 across a t-sweep ─────────────────────────
    @testset "L(t) ∈ [0, 1] sweep" begin
        m_0 = TFIM(; J=J, h=2.0)
        m_f = TFIM(; J=J, h=0.5)
        N = 24
        ts = range(0.0, 6.0; length=25)
        for t in ts
            L = QAtlas.fetch(m_f, LoschmidtEcho(:amplitude), OBC(N); initial=m_0, t=t)
            λ = QAtlas.fetch(m_f, LoschmidtRateFunction(), OBC(N); initial=m_0, t=t)
            @test -1e-10 ≤ L ≤ 1 + 1e-10
            @test λ ≥ -1e-10
        end
    end

    # ─── DQPT cusp at the analytical critical time (Infinite) ─────────────
    @testset "DQPT cusp: paramag → ferromag quench" begin
        h0, hf = 2.0, 0.5
        m_0 = TFIM(; J=J, h=h0)
        m_f = TFIM(; J=J, h=hf)

        kstar = _find_kstar(J, h0, hf)
        Λstar = _Λ(J, hf, kstar)
        t_c = π / (2 * Λstar)  # n = 0 critical time

        # The integrand has a log-divergence at k = k*, so λ(t_c) is finite
        # but the second derivative diverges — verify the cusp by checking
        # that λ at t_c sits in a local maximum relative to neighbours, and
        # that |slope_left − slope_right| at t_c is much larger than the
        # corresponding slope difference at a smooth control point.
        δ = 0.02
        λ_minus = QAtlas.fetch(
            m_f, LoschmidtRateFunction(), Infinite(); initial=m_0, t=t_c - δ
        )
        λ_at = QAtlas.fetch(m_f, LoschmidtRateFunction(), Infinite(); initial=m_0, t=t_c)
        λ_plus = QAtlas.fetch(
            m_f, LoschmidtRateFunction(), Infinite(); initial=m_0, t=t_c + δ
        )

        # Local maximum (cusp peak): λ(t_c) at least matches both neighbours.
        @test λ_at ≥ λ_minus - 1e-6
        @test λ_at ≥ λ_plus - 1e-6

        # Sign of slope flips across t_c.  Compare to a smooth region.
        t_smooth = 0.5 * t_c
        λ_sm_m = QAtlas.fetch(
            m_f, LoschmidtRateFunction(), Infinite(); initial=m_0, t=t_smooth - δ
        )
        λ_sm_p = QAtlas.fetch(
            m_f, LoschmidtRateFunction(), Infinite(); initial=m_0, t=t_smooth + δ
        )
        slope_smooth = (λ_sm_p - λ_sm_m) / (2δ)
        slope_left = (λ_at - λ_minus) / δ
        slope_right = (λ_plus - λ_at) / δ
        @test slope_left > 0
        @test slope_right < 0
        @test abs(slope_left - slope_right) > abs(slope_smooth) + 0.05
    end

    # ─── OBC large-N → Infinite agreement at off-DQPT t ───────────────────
    @testset "OBC N → ∞ matches Infinite (off-cusp)" begin
        h0, hf = 2.0, 0.5
        m_0 = TFIM(; J=J, h=h0)
        m_f = TFIM(; J=J, h=hf)
        # Pick t deliberately away from any t_n^* = π(n+1/2)/Λ_{k*}(h_f).
        kstar = _find_kstar(J, h0, hf)
        Λstar = _Λ(J, hf, kstar)
        t_c = π / (2 * Λstar)
        t = 0.5 * t_c  # well away from n=0,1 cusps

        λ_inf = QAtlas.fetch(m_f, LoschmidtRateFunction(), Infinite(); initial=m_0, t=t)
        λ_obc_small = QAtlas.fetch(m_f, LoschmidtRateFunction(), OBC(24); initial=m_0, t=t)
        λ_obc_large = QAtlas.fetch(m_f, LoschmidtRateFunction(), OBC(64); initial=m_0, t=t)

        @test λ_obc_small > 0
        @test λ_obc_large > 0
        # OBC bulk converges to Infinite as N grows.  After the factor-2
        # discretisation correction in `_tfim_loschmidt_obc_log_echo`, the
        # remaining Euler-Maclaurin error of the midpoint sum vs the
        # Infinite integral is O(1/N²) for smooth (off-cusp) integrands,
        # so N=64 should agree with the continuum limit to ≪ 1%.
        @test abs(λ_obc_large - λ_inf) < abs(λ_obc_small - λ_inf) + 0.05
        @test abs(λ_obc_large - λ_inf) < 1e-3
    end

    # ─── Pinned reference value: Infinite λ(t = 1.0; h_0=2, h_f=0.5) ──────
    @testset "Pinned λ(t=1, h_0=2, h_f=0.5)" begin
        m_0 = TFIM(; J=1.0, h=2.0)
        m_f = TFIM(; J=1.0, h=0.5)
        λ = QAtlas.fetch(m_f, LoschmidtRateFunction(), Infinite(); initial=m_0, t=1.0)
        @test 0.0 < λ < 1.0
        @test isfinite(λ)
        # Pinned value (atol 1e-8) — recorded after the first successful run
        # of QuadGK on the analytic integrand.  Update if the integrand
        # convention changes.
        @test λ ≈ 0.31693310885932685 atol = 1e-8
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "TFIM Loschmidt — verification cards" begin
    # No-quench identity: H_initial = H_final => |L(t)| = 1, rate λ(t) = 0
    # for every t (exact, independent of src).
    let m = TFIM(; J=1.0, h=1.5)
        for t in (0.5, 2.0, 7.3)
            verify(
                m,
                LoschmidtRateFunction(),
                Infinite();
                route=:limiting_case,
                fetch_kw=(; initial=m, t=t),
                independent=0.0,
                agree_within=1e-10,
                refs=["No-quench: H0 = Hf => λ(t) = 0 for all t"],
            )
        end
    end

    # t = 0: L(0) = ⟨ψ0|ψ0⟩ = 1 trivially for any quench pair
    verify(
        TFIM(; J=1.0, h=0.5),
        LoschmidtRateFunction(),
        Infinite();
        route=:limiting_case,
        fetch_kw=(; initial=TFIM(; J=1.0, h=2.0), t=0.0),
        independent=0.0,
        agree_within=1e-10,
        refs=["t=0: |L(0)| = 1 so the rate function λ(0) = 0"],
    )
end
