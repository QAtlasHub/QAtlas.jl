# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: TFIM fidelity susceptibility χ_F (issue #147)
#
# Validates the closed-form Bogoliubov-de Gennes implementation of
# `fetch(TFIM, FidelitySusceptibility(), bc)` for `bc ∈ {OBC, Infinite}`:
#
#  - Off-critical (h ≠ J): χ_F / L is finite and matches the closed-form
#    `1/(16(J² − h²))` (h < J) / `J²/(16 h² (h² − J²))` (h > J).
#  - Pinned numeric value at J = 1, h = 0.5 (Infinite): 1/12 (exact).
#  - L → ∞ extrapolation of OBC χ_F / N → Infinite value, off-critical.
#  - Critical scaling at h = J: χ_F ~ N² (i.e. χ_F / N ~ N — for the
#    1D-Ising universality class with d = 1, ν = 1).
#  - DomainError at h = J for `Infinite`.
#
# References:
#  - Gu, "Fidelity approach to quantum phase transitions",
#    Int. J. Mod. Phys. B 24, 4371 (2010).
#  - Damski, "Fidelity approach to quantum phase transitions in TFIM",
#    PRB 87, 165101 (2013).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM FidelitySusceptibility — Infinite, closed-form values" begin
    # Ordered phase (h < J): χ_F / L = 1 / (16 (J² − h²))
    @testset "ordered phase" begin
        for (J, h) in (
            (1.0, 0.0),
            (1.0, 0.1),
            (1.0, 0.5),
            (1.0, 0.9),
            (2.0, 0.5),
            (2.0, 1.5),
            (0.5, 0.2),
        )
            m = TFIM(J=J, h=h)
            χ = QAtlas.fetch(m, FidelitySusceptibility(), Infinite())
            expected = 1 / (16 * (J^2 - h^2))
            @test χ ≈ expected rtol = 1e-9
            @test isfinite(χ) && χ > 0
        end
    end

    # Disordered phase (h > J): χ_F / L = J² / (16 h² (h² − J²))
    @testset "disordered phase" begin
        for (J, h) in
            ((1.0, 1.5), (1.0, 2.0), (1.0, 5.0), (2.0, 3.0), (0.5, 1.0), (1.0, 1.01))
            m = TFIM(J=J, h=h)
            χ = QAtlas.fetch(m, FidelitySusceptibility(), Infinite())
            expected = J^2 / (16 * h^2 * (h^2 - J^2))
            @test χ ≈ expected rtol = 1e-9
            @test isfinite(χ) && χ > 0
        end
    end

    # Pinned value: J=1, h=0.5 → χ_F / L = 1 / (16 * 0.75) = 1/12
    @testset "pinned numeric value (J=1, h=0.5)" begin
        m = TFIM(J=1.0, h=0.5)
        χ = QAtlas.fetch(m, FidelitySusceptibility(), Infinite())
        @test χ ≈ 1 / 12 atol = 1e-12
        @test χ ≈ 0.08333333333333333 atol = 1e-12
    end

    # DomainError exactly at criticality
    @testset "DomainError at h = J" begin
        @test_throws DomainError QAtlas.fetch(
            TFIM(J=1.0, h=1.0), FidelitySusceptibility(), Infinite()
        )
        @test_throws DomainError QAtlas.fetch(
            TFIM(J=2.0, h=2.0), FidelitySusceptibility(), Infinite()
        )
        # ε-detuning OK
        χε = QAtlas.fetch(TFIM(J=1.0, h=1.0 - 1e-3), FidelitySusceptibility(), Infinite())
        @test isfinite(χε) && χε > 0
    end

    # Approach to criticality: divergence as h → J⁻
    @testset "approach to criticality (h → J⁻)" begin
        χ_h_close = QAtlas.fetch(TFIM(J=1.0, h=0.999), FidelitySusceptibility(), Infinite())
        χ_h_far = QAtlas.fetch(TFIM(J=1.0, h=0.5), FidelitySusceptibility(), Infinite())
        @test χ_h_close > 100 * χ_h_far
        # Linear divergence: χ_F / L ~ 1/(2 |h − J|) at leading order
        # (since χ_F / L = 1/(16 (J−h)(J+h)) ≈ 1/(32 (J − h)) as h → J).
        h_close = 0.999
        @test χ_h_close ≈ 1 / (16 * (1 - h_close^2)) rtol = 1e-9
    end
end

@testset "TFIM FidelitySusceptibility — OBC vs closed-form (off-critical)" begin
    # Sanity: OBC value finite, smooth, positive at every h ≠ J.
    m = TFIM(J=1.0, h=0.5)
    for N in (8, 16, 32, 64, 128)
        χ = QAtlas.fetch(m, FidelitySusceptibility(), OBC(N))
        @test isfinite(χ) && χ > 0
    end

    # L → ∞ extrapolation: χ_F(L) / L → Infinite value.
    # At h = 0.5 the gap is 2|J−h| = 1; finite-size corrections are O(1/N).
    @testset "convergence to thermodynamic limit (h=0.5)" begin
        target = QAtlas.fetch(m, FidelitySusceptibility(), Infinite())
        @test target ≈ 1 / 12 atol = 1e-12
        χ_per_site_64 = QAtlas.fetch(m, FidelitySusceptibility(), OBC(64)) / 64
        χ_per_site_128 = QAtlas.fetch(m, FidelitySusceptibility(), OBC(128)) / 128
        χ_per_site_256 = QAtlas.fetch(m, FidelitySusceptibility(), OBC(256)) / 256
        # Errors should shrink monotonically as N grows.
        err_64 = abs(χ_per_site_64 - target)
        err_128 = abs(χ_per_site_128 - target)
        err_256 = abs(χ_per_site_256 - target)
        @test err_128 < err_64
        @test err_256 < err_128
        # At N = 256 we should be within 5% of the thermodynamic value.
        @test err_256 < 0.05 * target
    end

    # Off-critical h = 0.3 also converges
    @testset "convergence at h=0.3" begin
        m2 = TFIM(J=1.0, h=0.3)
        target = QAtlas.fetch(m2, FidelitySusceptibility(), Infinite())
        χ_256 = QAtlas.fetch(m2, FidelitySusceptibility(), OBC(256)) / 256
        @test abs(χ_256 - target) < 0.05 * target
    end

    # Disordered side: h = 2.0
    @testset "convergence at h=2.0 (disordered)" begin
        m3 = TFIM(J=1.0, h=2.0)
        target = QAtlas.fetch(m3, FidelitySusceptibility(), Infinite())
        @test target ≈ 1 / 192 atol = 1e-12
        χ_256 = QAtlas.fetch(m3, FidelitySusceptibility(), OBC(256)) / 256
        @test abs(χ_256 - target) < 0.05 * target
    end
end

@testset "TFIM FidelitySusceptibility — OBC critical scaling" begin
    # At h = J the 1D Ising universality class predicts χ_F / L ~ |h − h_c|^{−1}
    # in the thermodynamic limit; finite-size scaling at h = h_c gives
    # χ_F(N) ~ N² (i.e. χ_F / N ~ N).  We test that χ_F(N) / N² approaches
    # a constant as N grows (relative variation below ~2% across N = 64,
    # 128, 256).  This is the standard ν-scaling consistency check
    # χ_F ~ N^{2/ν − d + 2} with ν = 1, d = 1 (Damski 2013, Sec. III.B).
    m = TFIM(J=1.0, h=1.0)
    ratios = Float64[]
    Ns = (64, 128, 256)
    for N in Ns
        χ = QAtlas.fetch(m, FidelitySusceptibility(), OBC(N))
        push!(ratios, χ / N^2)
    end
    # Successive ratios stabilise (asymptotic constant).
    @test isapprox(ratios[2], ratios[1]; rtol=0.02)
    @test isapprox(ratios[3], ratios[2]; rtol=0.01)
    # All ratios positive and around the empirical universal value
    # (~ 0.0177–0.0179 by direct numerical measurement at finite N).
    for r in ratios
        @test 0.005 < r < 0.05
    end
end

@testset "TFIM FidelitySusceptibility — internal consistency" begin
    # χ_F is invariant under h → -h for any J > 0 (the model has σˣ → -σˣ
    # parity; see core Hamiltonian).  We test the Infinite case.
    for h in (0.3, 0.7, 1.5, 3.0)
        m_plus = TFIM(J=1.0, h=h)
        m_minus = TFIM(J=1.0, h=(-h))
        @test QAtlas.fetch(m_plus, FidelitySusceptibility(), Infinite()) ≈
            QAtlas.fetch(m_minus, FidelitySusceptibility(), Infinite()) atol = 1e-10
    end

    # Per-site option: per_site=true returns χ_F / N.
    m = TFIM(J=1.0, h=0.5)
    χ_total = QAtlas.fetch(m, FidelitySusceptibility(), OBC(32))
    χ_persite = QAtlas.fetch(m, FidelitySusceptibility(), OBC(32); per_site=true)
    @test χ_persite ≈ χ_total / 32 atol = 1e-12
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "TFIM FidelitySusceptibility — verification cards" begin
    # Bogoliubov-de Gennes closed form for the per-site fidelity
    # susceptibility (independent re-derivation, not src):
    #   ordered (h < J):    χ_F / L = 1 / (16 (J² − h²))
    #   disordered (h > J): χ_F / L = J² / (16 h² (h² − J²))
    let J = 1.0, h = 0.5
        verify(
            TFIM(; J=J, h=h),
            FidelitySusceptibility(),
            Infinite();
            route=:second_closed_form,
            independent=1 / (16 * (J^2 - h^2)),
            agree_within=1e-9,
            refs=["BdG closed form: χ_F/L = 1/(16(J²−h²)) ordered phase (= 1/12 at J=1,h=1/2)"],
        )
    end
    let J = 1.0, h = 2.0
        verify(
            TFIM(; J=J, h=h),
            FidelitySusceptibility(),
            Infinite();
            route=:second_closed_form,
            independent=J^2 / (16 * h^2 * (h^2 - J^2)),
            agree_within=1e-9,
            refs=["BdG closed form: χ_F/L = J²/(16h²(h²−J²)) disordered phase"],
        )
    end
end
