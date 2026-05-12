# =============================================================================
# Tests for the TFIM real-time dynamics module.
#
# Three layers of validation:
#   1. Equal-time and unequal-time correlators against small-N exact
#      diagonalization (ED) — pure correctness check.
#   2. Critical-point spatial decay  ⟨σᶻ_0 σᶻ_r⟩ ~ r^{-1/4}
#      and temporal decay  ⟨σᶻ(t) σᶻ(0)⟩ ~ t^{-1/4}
#      from the Ising CFT (scaling dimension Δ_σ = 1/8).
#   3. Off-critical (disordered phase, h > J) exponential spatial decay
#      with the analytically known correlation length.
# =============================================================================

# Average of ⟨σᶻ_i σᶻ_{i+r}⟩ over a band of bulk r values, used as a
# numerical estimate of the long-range plateau (in the ordered phase).
function mean_far(N::Int, J::Float64, h::Float64, i0::Int)
    rs = filter(r -> i0 + r ≤ N - 5, [10, 12, 14, 16])
    isempty(rs) && error("no bulk r values fit; choose larger N or smaller i0")
    s = 0.0
    for r in rs
        s += real(
            QAtlas.fetch(
                TFIM(; J=J, h=h),
                ZZCorrelation{:dynamic}(),
                OBC(; N=N);
                i=i0,
                j=i0 + r,
                t=0.0,
            ),
        )
    end
    return s / length(rs)
end

# ============================================================================

@testset "TFIM dynamics" begin

    # ---- Layer 1: ED comparison, small N ------------------------------------
    @testset "ED comparison N=4, J=1, h=0.5" begin
        N, J, h = 4, 1.0, 0.5
        H = _build_tfim_dense(N, J, h)
        E, V = eigen(H)
        gs = V[:, 1]

        # Sanity: GS energy matches QAtlas BdG.
        @test E[1] ≈ QAtlas.fetch(TFIM(; J=J, h=h), Energy(), OBC(; N=N)) atol=1e-10

        # σz σz at t = 0
        for i in 1:N, j in i:N
            ED = real(gs' * _op_site(_SZ, i, N) * _op_site(_SZ, j, N) * gs)
            QAT = real(
                QAtlas.fetch(
                    TFIM(; J=J, h=h), ZZCorrelation{:dynamic}(), OBC(; N=N); i=i, j=j, t=0.0
                ),
            )
            @test QAT ≈ ED atol=1e-10
        end

        # σx σx at t = 0
        for i in 1:N, j in i:N
            ED = real(gs' * _op_site(_SX, i, N) * _op_site(_SX, j, N) * gs)
            QAT = real(
                QAtlas.fetch(
                    TFIM(; J=J, h=h), XXCorrelation{:dynamic}(), OBC(; N=N); i=i, j=j, t=0.0
                ),
            )
            @test QAT ≈ ED atol=1e-10
        end

        # Unequal-time σz σz, several (i, j, t) — full complex check
        Udag = V'  # = V^†
        for (i, j, t) in [(1, 1, 0.3), (2, 3, 0.7), (1, 4, 1.5), (2, 4, 2.0)]
            Ut = V * Diagonal(exp.(-im * E * t)) * Udag
            ED =
                exp(im * E[1] * t) *
                (gs' * _op_site(_SZ, i, N) * Ut * _op_site(_SZ, j, N) * gs)
            QAT = QAtlas.fetch(
                TFIM(; J=J, h=h), ZZCorrelation{:dynamic}(), OBC(; N=N); i=i, j=j, t=t
            )
            @test real(QAT) ≈ real(ED) atol=1e-10
            @test imag(QAT) ≈ imag(ED) atol=1e-10
        end
    end

    # =====================================================================
    # Layer 2: rigorous analytical results vs. finite-size / finite-r data.
    # The QAtlas correlator is *exact* for finite N, so the test pattern is:
    #
    #   ‖ exact(QAtlas, N, r) − analytic_TL ‖ → 0  as  N → ∞ or r → ∞
    #
    # We don't fit phenomenological tolerances; we check that the deviation
    # *shrinks monotonically* as the scaling parameter is increased toward
    # the thermodynamic limit (or large-r asymptotic).
    # =====================================================================

    # ---- Layer 2a: Pfeuty exact magnetisation in the ordered phase ---------
    # Pfeuty (1970) gives the exact thermodynamic-limit magnetisation
    #   M(h) = (1 - (h/J)^2)^{1/8}     (h ≤ J)
    # so the long-range plateau of ⟨σᶻ_i σᶻ_j⟩ in the ordered phase equals
    #   M² = (1 - (h/J)²)^{1/4}.
    # On a finite OBC chain the GS is Z₂-symmetric (⟨σᶻ⟩ = 0), so the
    # connected and disconnected correlators coincide and the bulk-bulk
    # plateau converges to M² *very* fast (corrections are exponential in
    # N times the inverse-gap scale).  We test the convergence on a
    # sequence of N values.
    @testset "ordered phase: Pfeuty M² convergence" begin
        J = 1.0
        # Use h values where the convergence is visible above machine
        # precision (closer to criticality ⇒ slower convergence).
        for h in (0.6, 0.7, 0.8)
            M2_exact = (1 - (h / J)^2)^(1 / 4)
            errs = Float64[]
            for N in (30, 50, 80)
                i0 = N ÷ 4
                far = mean_far(N, J, h, i0)
                push!(errs, abs(far - M2_exact))
            end
            # Convergence: largest-N error must shrink by at least 10× over
            # the N range, and the largest-N error must be < 1e-4.  Both
            # are very loose for the actual exponential convergence.
            @test errs[end] < errs[1] / 10
            @test errs[end] < 1e-4
        end
    end

    # ---- Layer 2b: Critical CFT exponent — finite-N scaling toward -1/4 ----
    # Ising CFT predicts ⟨σ_0 σ_r⟩ ~ r^{-1/4} in the TL (Δ_σ = 1/8).  At
    # finite N with OBC the *effective* slope estimated from a doubling
    # ratio (r, 2r) at fixed bulk centre converges to -1/4 as N grows.  We
    # require monotonic decrease of the deviation and that it falls below
    # 0.10 by the largest N tested.
    @testset "criticality: scaling toward CFT exponent -1/4" begin
        J = 1.0
        h = 1.0
        Ns = (80, 160, 240)
        slope_exact = -1 / 4
        errs = Float64[]
        for N in Ns
            i0 = N ÷ 4
            v1 = real(
                QAtlas.fetch(
                    TFIM(; J=J, h=h),
                    ZZCorrelation{:dynamic}(),
                    OBC(; N=N);
                    i=i0,
                    j=i0 + 10,
                    t=0.0,
                ),
            )
            v2 = real(
                QAtlas.fetch(
                    TFIM(; J=J, h=h),
                    ZZCorrelation{:dynamic}(),
                    OBC(; N=N);
                    i=i0,
                    j=i0 + 20,
                    t=0.0,
                ),
            )
            slope = log2(v2 / v1)         # since 20/10 = 2
            push!(errs, abs(slope - slope_exact))
        end
        @test issorted(errs; rev=true)    # monotone convergence
        @test errs[end] < 0.10            # within 0.10 of the CFT value at N=240
    end

    # ---- Layer 2c: Disordered exact correlation length — r-scaling ---------
    # In the disordered phase (h > J) the *exact* thermodynamic-limit inverse
    # correlation length follows from the BdG dispersion
    #   ω(k) = 2 √(J² + h² - 2 J h cos k)
    # via continuation k → i/ξ:
    #   ξ⁻¹ = arccosh((h² + J²) / (2 h J)).
    # The local log-derivative of ⟨σᶻ_0 σᶻ_r⟩ approaches -ξ⁻¹ from below as
    # r grows; we verify monotone convergence on a sequence of (r, r+Δ)
    # pairs.  Floating-point precision sets an upper limit on r (the
    # correlator decays exponentially), so we stop at r = 24 where
    # |C| ~ 10⁻⁸ for our test parameters.
    @testset "disordered: r-scaling toward exact ξ⁻¹" begin
        J = 1.0
        for h in (1.5, 2.0)
            slope_exact = -acosh((h^2 + J^2) / (2 * h * J))
            N = 200
            i0 = 50
            pairs = [(8, 10), (12, 14), (16, 20), (20, 24)]
            errs = Float64[]
            for (r1, r2) in pairs
                v1 = real(
                    QAtlas.fetch(
                        TFIM(; J=J, h=h),
                        ZZCorrelation{:dynamic}(),
                        OBC(; N=N);
                        i=i0,
                        j=i0 + r1,
                        t=0.0,
                    ),
                )
                v2 = real(
                    QAtlas.fetch(
                        TFIM(; J=J, h=h),
                        ZZCorrelation{:dynamic}(),
                        OBC(; N=N);
                        i=i0,
                        j=i0 + r2,
                        t=0.0,
                    ),
                )
                slope = (log(abs(v2)) - log(abs(v1))) / (r2 - r1)
                push!(errs, abs(slope - slope_exact))
            end
            @test issorted(errs; rev=true)   # monotone convergence
            @test errs[end] < 0.025          # exact within 2.5 % at r ≈ 22
        end
    end

    # ---- Layer 2d: criticality — temporal envelope decay -------------------
    # At fixed bulk site, |⟨σᶻ(t) σᶻ(0)⟩|_GS decays slowly at criticality.
    # Inside the boundary-bounce time it should be monotone (after the
    # first sign of the lattice short-time wiggle has passed) and stay
    # well above any exponential-gap bound (the gap closes).
    @testset "criticality: temporal envelope decay" begin
        N, J, h = 60, 1.0, 1.0
        i0 = N ÷ 2
        ts = [2.0, 4.0, 8.0]
        envelope = Float64[
            abs(
                QAtlas.fetch(
                    TFIM(; J=J, h=h), ZZCorrelation{:dynamic}(), OBC(; N=N); i=i0, j=i0, t=t
                ),
            ) for t in ts
        ]
        @test issorted(envelope; rev=true)
        @test envelope[end] > 0.05        # not exponential
    end

    # ---- Sanity: equal-position autocorrelator at t=0 is 1 ------------------
    @testset "(σᶻ_i)² = 1 sanity" begin
        N, J, h = 12, 1.0, 0.7
        for i in (1, 4, 8, 12)
            v = QAtlas.fetch(
                TFIM(; J=J, h=h), ZZCorrelation{:dynamic}(), OBC(; N=N); i=i, j=i, t=0.0
            )
            @test real(v) ≈ 1.0 atol=1e-10
            @test abs(imag(v)) < 1e-10
        end
    end

    # ---- spreading convenience function ------------------------------------
    @testset "sz_sz_spreading shape and consistency" begin
        N, J, h = 10, 1.0, 0.5
        center = 5
        times = [0.0, 0.5, 1.0]
        C = QAtlas.fetch(
            TFIM(; J=J, h=h),
            ZZCorrelation{:lightcone}(),
            OBC(; N=N);
            center=center,
            times=times,
        )
        @test size(C) == (3, N)

        # Each entry should match an individual call.
        for (it, t) in enumerate(times), ix in 1:N
            ref = QAtlas.fetch(
                TFIM(; J=J, h=h), ZZCorrelation{:dynamic}(), OBC(; N=N); i=ix, j=center, t=t
            )
            @test C[it, ix] ≈ ref atol=1e-10
        end
    end
end
