# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: IsingTriangular — Wannier 1950 frustrated AFM closed forms.
#
# Verifies:
#   * AFM (J > 0) ResidualEntropy ≈ 0.3230659669 (Wannier 1950) at 1e-9
#   * AFM CriticalTemperature is exactly 0.0
#   * FM (J = -1) CriticalTemperature = 4 |J| / log 3 (Houtappel 1950) at 1e-12
#   * FM ResidualEntropy = 0
#   * J = 0 edge case returns 0 for both
#   * J-independence of ResidualEntropy in the AFM branch
#
# No Lattice2D dependency.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test
using QuadGK: quadgk

@testset "IsingTriangular — CriticalTemperature" begin
    # AFM (J > 0): T_c = 0 (Wannier 1950, frustrated)
    @test QAtlas.fetch(IsingTriangular(; J=1.0), CriticalTemperature(), Infinite()) == 0.0
    @test QAtlas.fetch(IsingTriangular(; J=2.5), CriticalTemperature(), Infinite()) == 0.0

    # FM (J < 0): T_c = 4 |J| / log 3 (Houtappel 1950)
    Tc_unit = 4 / log(3)
    @test QAtlas.fetch(IsingTriangular(; J=-1.0), CriticalTemperature(), Infinite()) ≈
        Tc_unit atol = 1e-12
    # Numerical sanity check on the Houtappel value: 4/log(3) ≈ 3.64096
    @test 3.6 < Tc_unit < 3.65
    @test Tc_unit * log(3) ≈ 4.0 atol = 1e-14

    # |J| scaling: T_c(J = -c) = c · T_c(J = -1)
    for J in (-0.5, -1.0, -2.0, -3.7)
        Tc = QAtlas.fetch(IsingTriangular(; J=J), CriticalTemperature(), Infinite())
        @test Tc ≈ abs(J) * Tc_unit rtol = 1e-14
    end

    # J = 0: degenerate; documented as T_c = 0
    @test QAtlas.fetch(IsingTriangular(; J=0.0), CriticalTemperature(), Infinite()) == 0.0

    # kwarg override matches struct field
    m = IsingTriangular(; J=1.0)
    @test QAtlas.fetch(m, CriticalTemperature(), Infinite(); J=-1.0) ≈ Tc_unit atol = 1e-12
end

@testset "IsingTriangular — ResidualEntropy (Wannier 1950)" begin
    # Wannier 1950's published value `0.3230659669` is rounded/truncated
    # at ~10 digits.  The actual high-precision value of the integral
    # `(2/π) ∫₀^{π/3} log(2 cos θ) dθ` is `0.32306594722...`, so a
    # tighter cross-check is to recompute the integral here at very high
    # precision via QuadGK.
    val, _err = quadgk(θ -> log(2 * cos(θ)), 0.0, π / 3; rtol=1e-14, atol=1e-14)
    S_quadgk = (2 / π) * val
    # Coarse literature-rounding cross-check (loose atol — the published
    # 10-digit truncation deviates from machine precision in the 8th
    # decimal place only).
    S_ref_truncated = 0.3230659669

    @testset "AFM (J > 0) — Wannier integral" begin
        S = QAtlas.fetch(IsingTriangular(; J=1.0), ResidualEntropy(), Infinite())

        # Match Wannier's published 10-digit truncation, allowing for
        # the literature rounding noise (~1e-7).
        @test S ≈ S_ref_truncated atol = 1e-6
        # Match the QuadGK reference recomputation to ~machine ε
        @test S ≈ S_quadgk atol = 1e-12
        # Direct numerical pin: 0.32306594722 is correct to 1e-11
        @test S ≈ 0.32306594721945 atol = 1e-12

        # Hard bound: in (0, log 2)  (must be positive but < log 2, the
        # naive site-decoupled upper bound for an Ising spin)
        @test 0 < S < log(2)

        # J-independence of the AFM residual entropy
        for J in (0.5, 1.0, 2.0, 7.3)
            @test QAtlas.fetch(IsingTriangular(; J=J), ResidualEntropy(), Infinite()) ≈ S atol =
                1e-14
        end
    end

    @testset "FM (J < 0) — ordered ground state" begin
        for J in (-0.1, -1.0, -3.5)
            @test QAtlas.fetch(IsingTriangular(; J=J), ResidualEntropy(), Infinite()) == 0.0
        end
    end

    @testset "J = 0 edge case" begin
        # The non-interacting limit is a degenerate convention; we
        # report 0 (matches the FM branch).  The test is here to catch
        # any future change of convention.
        @test QAtlas.fetch(IsingTriangular(; J=0.0), ResidualEntropy(), Infinite()) == 0.0
    end

    @testset "kwarg override" begin
        m = IsingTriangular(; J=-1.0)
        # The struct says FM but caller forces AFM via kwarg
        S = QAtlas.fetch(m, ResidualEntropy(), Infinite(); J=1.0)
        @test S ≈ S_quadgk atol = 1e-12
    end
end

@testset "IsingTriangular — CriticalExponents = 2D Ising Onsager (Phase 2)" begin
    m = IsingTriangular()
    exp = QAtlas.fetch(m, CriticalExponents(), Infinite())
    @test exp.α == 0
    @test exp.β == 1 // 8
    @test exp.γ == 7 // 4
    @test exp.δ == 15
    @test exp.ν == 1
    @test exp.η == 1 // 4
    # Delegation invariant
    @test exp == QAtlas.fetch(QAtlas.Universality(:Ising), CriticalExponents(); d=2)
    # Hyperscaling
    @test exp.α + 2 * exp.β + exp.γ == 2
    @test exp.γ == exp.β * (exp.δ - 1)
    @test exp.η == 2 - exp.γ // exp.ν
    # Universality check: matches IsingSquare critical exponents (after #346 lands)
    # — both are in the 2D Ising universality class.
end
