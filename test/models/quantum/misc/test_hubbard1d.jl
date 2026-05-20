# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: Hubbard1D Lieb–Wu Bethe ansatz half-filling closed form.
#
# Run via:
#
#   julia --project=test test/standalone/test_hubbard1d.jl
#
# Cases (issue #153 phase 1):
#   1. Pinned analytical E₀/N at U/t = 4   (numerical reference value).
#   2. E₀/N → -4t/π                as U/t → 0.
#   3. E₀/N → -4 t² log 2 / U      as U/t → ∞.
#   4. Charge gap > 0 at any U > 0 (Mott — no transition in 1D Hubbard).
#   5. Charge gap → 0              as U → 0.
#   6. Charge gap → U - 4t         as U → ∞.
#   7. Spin gap = 0                exactly.
#   8. Off half filling (μ = U/3)  raises DomainError.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

const _PINNED_E0_U4 = -0.5737293678984494  # quadgk rtol=1e-12 at t=1, U=4
const _PINNED_DC_U4 = 1.2867270220129354   # quadgk rtol=1e-12 at t=1, U=4

@testset "Hubbard1D — Lieb–Wu half-filling Phase 1" begin
    @testset "Constructor + half-filling default" begin
        m = Hubbard1D()
        @test m.t == 1.0
        @test m.U == 4.0
        @test m.μ == 2.0   # half filling (U/2)

        m2 = Hubbard1D(; t=2.0, U=8.0, μ=4.0)
        @test (m2.t, m2.U, m2.μ) == (2.0, 8.0, 4.0)
    end

    @testset "Pinned E₀/N at U/t = 4" begin
        m = Hubbard1D(; t=1.0, U=4.0, μ=2.0)
        e0 = QAtlas.fetch(m, GroundStateEnergyDensity(), Infinite())
        @test e0 ≈ _PINNED_E0_U4 atol = 1e-8
    end

    @testset "E₀/N → -4t/π as U/t → 0 (free 1D fermion)" begin
        # At U/t = 0.01 the Lieb–Wu integral should be very close to
        # the free-fermion value -4t/π.
        m = Hubbard1D(; t=1.0, U=0.01, μ=0.005)
        e0 = QAtlas.fetch(m, GroundStateEnergyDensity(), Infinite())
        @test e0 ≈ -4.0 / π atol = 1e-2
    end

    @testset "E₀/N → -4 t² log(2) / U as U/t → ∞ (Heisenberg AFM)" begin
        # At U/t = 100 the half-filled Hubbard chain reduces to the
        # Heisenberg AFM with effective coupling J = 4 t²/U; the GS
        # energy density is -J log 2 = -4 t² log 2 / U.
        m = Hubbard1D(; t=1.0, U=100.0, μ=50.0)
        e0 = QAtlas.fetch(m, GroundStateEnergyDensity(), Infinite())
        e0_asymp = -4.0 * 1.0^2 * log(2.0) / 100.0
        @test isapprox(e0, e0_asymp; rtol=0.05)
    end

    @testset "Charge gap pinned at U/t = 4" begin
        m = Hubbard1D(; t=1.0, U=4.0, μ=2.0)
        Δc = QAtlas.fetch(m, ChargeGap(), Infinite())
        @test Δc ≈ _PINNED_DC_U4 atol = 1e-8
        @test Δc > 0  # Mott
    end

    @testset "Charge gap > 0 for any U > 0 (Mott — no 1D transition)" begin
        for U in (0.5, 1.0, 4.0, 16.0)
            m = Hubbard1D(; t=1.0, U=U, μ=U / 2)
            Δc = QAtlas.fetch(m, ChargeGap(), Infinite())
            @test Δc > 0
        end
    end

    @testset "Charge gap → 0 as U → 0" begin
        # At U = 0.1 the Mott gap is exponentially small (~ exp(-2π t/U)).
        m = Hubbard1D(; t=1.0, U=0.1, μ=0.05)
        Δc = QAtlas.fetch(m, ChargeGap(), Infinite())
        @test Δc < 0.05
        @test isapprox(Δc, 0.0; atol=0.05)
    end

    @testset "Charge gap → U - 4t as U → ∞" begin
        # At U = 100, t = 1, Δ_c → 100 - 4 + 8 log 2 / 100 ≈ 96.055.
        # Atol 0.1 is comfortable.
        m = Hubbard1D(; t=1.0, U=100.0, μ=50.0)
        Δc = QAtlas.fetch(m, ChargeGap(), Infinite())
        @test isapprox(Δc, 100.0 - 4.0; atol=0.1)
    end

    @testset "Spin gap = 0 exactly" begin
        for (t, U) in ((1.0, 4.0), (1.0, 0.5), (2.0, 8.0), (1.0, 100.0))
            m = Hubbard1D(; t=t, U=U, μ=U / 2)
            Δs = QAtlas.fetch(m, SpinGap(), Infinite())
            @test Δs == 0.0
        end
    end

    @testset "Off half-filling raises DomainError" begin
        # μ = U/3 (≠ U/2) — Phase 2 territory.
        m = Hubbard1D(; t=1.0, U=6.0, μ=2.0)  # U/2 = 3.0, μ = 2.0
        @test_throws DomainError QAtlas.fetch(m, GroundStateEnergyDensity(), Infinite())
        @test_throws DomainError QAtlas.fetch(m, ChargeGap(), Infinite())
        @test_throws DomainError QAtlas.fetch(m, SpinGap(), Infinite())
    end
end

@testset "Hubbard1D — LuttingerParameter free-fermion U=0 (Phase 2)" begin
    # Default Hubbard1D has U=4 — construct explicit U=0
    K = QAtlas.fetch(Hubbard1D(; U=0.0), LuttingerParameter(), Infinite())
    @test K == 1.0
    # Various t values — K is t-independent (universal LL value)
    for t in (0.5, 1.0, 3.0)
        @test QAtlas.fetch(Hubbard1D(; t=t, U=0.0), LuttingerParameter(), Infinite()) == 1.0
    end
end

@testset "Hubbard1D — LuttingerParameter U ≠ 0 throws DomainError (Phase 2 deferral)" begin
    @test_throws DomainError QAtlas.fetch(
        Hubbard1D(; U=4.0), LuttingerParameter(), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        Hubbard1D(; U=-1.0), LuttingerParameter(), Infinite()
    )
    # Regression: strict iszero(U) — tiny non-zero U must NOT silently return K=1
    # (K(U) is not analytic-equal to 1 for any U ≠ 0; Lieb-Wu integrals required).
    @test_throws DomainError QAtlas.fetch(
        Hubbard1D(; U=1e-13), LuttingerParameter(), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        Hubbard1D(; U=-1e-13), LuttingerParameter(), Infinite()
    )
end

# ── additional verification cards (#381 batch) ─────────────────────────────
# Scope note: GroundStateEnergyDensity at U → 0 was deliberately omitted
# from this batch — the src `_hubbard1d_e0` returns -4t²/π (not the
# textbook free-fermion value -4t/π) at the U → 0 half-filling limit;
# the prefactor discrepancy is a separate src issue from corroboration.
@testset "Hubbard1D — additional Lieb-Wu cards (#381 batch)" begin
    # SpinGap/Infinite: rigorous Lieb-Wu (1968) — spinon sector is gapless
    # at any U > 0 by SU(2) symmetry + Bethe ansatz. Δ_s = 0 exactly.
    for (t, U) in ((1.0, 0.5), (1.0, 4.0), (2.0, 8.0))
        verify(
            Hubbard1D(; t=t, U=U, μ=U/2),
            SpinGap(),
            Infinite();
            route=:second_closed_form,
            independent=0.0,
            agree_within=0,
            refs=["Lieb-Wu 1968 (Bethe ansatz): Δ_s = 0 for all U > 0 (gapless spinons by SU(2) symmetry)"],
        )
    end

    # ChargeGap/Infinite, U → 0 limit: the Mott gap is exponentially small,
    # Δ_c ∝ exp(-2π t / U). Using U = 0.3 t gives Δ_c ≈ exp(-2π/0.3) ≈ 5e-10,
    # which IS representable in Float64 (the earlier U = 0.05 t put the gap
    # at ~10⁻⁵⁵, deep under Float64 normal range — that would test
    # underflow-to-zero, not the limiting-case physics). The
    # route=:limiting_case marker (not :second_closed_form like the
    # sibling cards) signals that Δ_c → 0 is an asymptotic statement,
    # not an exact closed form at finite U.
    for t in (0.5, 1.0, 2.0)
        U = 3e-1 * t
        verify(
            Hubbard1D(; t=t, U=U, μ=U/2),
            ChargeGap(),
            Infinite();
            route=:limiting_case,
            independent=0.0,
            agree_within=1e-4,
            refs=["Lieb-Wu 1968: Δ_c → 0 as U → 0 with exponential form Δ_c ∝ exp(-2π t / U)"],
        )
    end
end

