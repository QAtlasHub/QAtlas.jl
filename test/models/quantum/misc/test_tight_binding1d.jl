# ─────────────────────────────────────────────────────────────────────────────
# TightBinding1D — 1D non-interacting spinless-fermion chain (Phase 1).
#
# Targeted run (skips Pkg.test()):
#   julia --project=test test/models/quantum/misc/test_tight_binding1d.jl
#
# Coverage:
#   • Energy{:per_site} at Infinite — half-filling -2/π, t-linearity,
#     band-edge limits μ = ±2t, deep empty / full band, generic μ = 1
#     (k_F = 2π/3 → E/N = -√3/π - 2/3).
#   • MassGap at Infinite — gapless interior, exact closed form on the
#     insulating side, t-scaling check.
#   • FermiVelocity at Infinite — v_F = 2t at half-filling, √3 at μ = 1,
#     t-linearity, returns 0 in the gapped phase (|μ| ≥ 2t).
#   • Constructor — rejects t ≤ 0 with DomainError.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test
using QAtlas: TightBinding1D, Energy, MassGap, FermiVelocity, NMRSpinRelaxationRate, Infinite, fetch

@testset "TightBinding1D" begin

    # ──────────────────────────── Energy{:per_site} ────────────────────────────
    @testset "Energy{:per_site} — closed-form ground-state energy density" begin
        # Half-filling, t = 1, μ = 0: k_F = π/2 → E/N = -(2/π) sin(π/2) - 0 = -2/π.
        @test QAtlas.fetch(TightBinding1D(), Energy{:per_site}(), Infinite()) ≈ -2 / π

        # Linear in t at fixed μ/t: t = 3, μ = 0 → E/N = -6/π.
        @test QAtlas.fetch(TightBinding1D(; t=3.0), Energy{:per_site}(), Infinite()) ≈
            -6 / π

        # Band-edge limits.  At μ = -2t the partial-filling integral collapses
        # smoothly to 0 (empty band); at μ = +2t to -μ (full band).
        @test QAtlas.fetch(
            TightBinding1D(; t=1.0, μ=-2.0), Energy{:per_site}(), Infinite()
        ) == 0.0
        @test QAtlas.fetch(
            TightBinding1D(; t=1.0, μ=2.0), Energy{:per_site}(), Infinite()
        ) ≈ -2.0

        # Deep insulating regimes (μ outside ±2t): empty / full band branches.
        @test QAtlas.fetch(
            TightBinding1D(; t=1.0, μ=3.0), Energy{:per_site}(), Infinite()
        ) == -3.0                                                # full band
        @test QAtlas.fetch(
            TightBinding1D(; t=1.0, μ=-5.0), Energy{:per_site}(), Infinite()
        ) == 0.0                                                 # empty band

        # Generic μ = 1, t = 1 — k_F = arccos(-1/2) = 2π/3:
        #   E/N = -(2/π) sin(2π/3) - (1/π)(2π/3) = -√3/π - 2/3 ≈ -1.21790...
        @test QAtlas.fetch(
            TightBinding1D(; t=1.0, μ=1.0), Energy{:per_site}(), Infinite()
        ) ≈ -sqrt(3) / π - 2 / 3
    end

    # ──────────────────────────────── MassGap ──────────────────────────────────
    @testset "MassGap — Δ = max(0, |μ| - 2t)" begin
        # Half-filling and shallow metal — gapless.
        @test QAtlas.fetch(TightBinding1D(; t=1.0, μ=0.0), MassGap(), Infinite()) == 0.0
        @test QAtlas.fetch(TightBinding1D(; t=1.0, μ=1.5), MassGap(), Infinite()) == 0.0

        # Lifshitz transition — band edge exactly at chemical potential.
        @test QAtlas.fetch(TightBinding1D(; t=1.0, μ=2.0), MassGap(), Infinite()) == 0.0

        # Insulating: μ = ±3 with t = 1 → gap = |μ| - 2.
        @test QAtlas.fetch(TightBinding1D(; t=1.0, μ=3.0), MassGap(), Infinite()) == 1.0
        @test QAtlas.fetch(TightBinding1D(; t=1.0, μ=-5.0), MassGap(), Infinite()) == 3.0

        # Bandwidth scaling: at t = 2 the metallic window is |μ| < 4,
        # so μ = 5 gives gap 5 - 4 = 1.
        @test QAtlas.fetch(TightBinding1D(; t=2.0, μ=5.0), MassGap(), Infinite()) == 1.0
    end

    # ─────────────────────────────── FermiVelocity ─────────────────────────────
    @testset "FermiVelocity — v_F = 2t √(1 - μ²/(4t²))" begin
        # Half-filling: k_F = π/2, sin(k_F) = 1 → v_F = 2t.
        @test QAtlas.fetch(TightBinding1D(), FermiVelocity(), Infinite()) ≈ 2.0
        @test QAtlas.fetch(TightBinding1D(; t=3.0), FermiVelocity(), Infinite()) ≈ 6.0

        # μ = 1, t = 1: sin(2π/3) = √3/2 → v_F = √3.
        @test QAtlas.fetch(TightBinding1D(; t=1.0, μ=1.0), FermiVelocity(), Infinite()) ≈
            sqrt(3)

        # No Fermi surface in the insulating regime — both signs of μ
        # return 0 by convention (gapped phase).
        @test QAtlas.fetch(TightBinding1D(; t=1.0, μ=2.5), FermiVelocity(), Infinite()) ==
            0.0
        @test QAtlas.fetch(TightBinding1D(; t=1.0, μ=-3.0), FermiVelocity(), Infinite()) ==
            0.0
    end

    # ────────────────────────── Constructor validation ─────────────────────────
    @testset "Constructor rejects t ≤ 0" begin
        @test_throws DomainError TightBinding1D(; t=0.0)
        @test_throws DomainError TightBinding1D(; t=-1.0)
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "TightBinding1D — verification cards" begin
    verify(
        TightBinding1D(; t=1.0, μ=0.0),
        Energy(:per_site),
        Infinite();
        route=:second_closed_form,
        independent=-2 / pi,
        agree_within=1e-9,
        refs=["Half-filled tight-binding chain: e0 = -2t/pi"],
    )
    verify(
        TightBinding1D(; t=1.0, μ=0.0),
        MassGap(),
        Infinite();
        route=:second_closed_form,
        independent=0.0,
        agree_within=1e-10,
        refs=["Half filling: gapless Fermi surface"],
    )
    verify(
        TightBinding1D(; t=1.0, μ=3.0),
        MassGap(),
        Infinite();
        route=:second_closed_form,
        independent=1.0,
        agree_within=1e-9,
        refs=["Band insulator mu > 2t: gap = |mu| - 2t"],
    )
end
# ── additional verification cards (#381 batch) ─────────────────────────────
@testset "TightBinding1D — additional verification cards (#381 batch)" begin
    # Free-fermion v_F = 2 t sin(k_F); default t=1, μ=0 ⇒ k_F = π/2 ⇒ v_F = 2.
    verify(
        TightBinding1D(),
        FermiVelocity(),
        Infinite();
        route=:second_closed_form,
        independent=2.0,
        agree_within=1e-12,
        refs=["Ashcroft-Mermin 1976: v_F = 2 t sin(k_F); μ=0 half-filling"],
    )

    # NMR spin relaxation rate low-T Korringa limit with regularization eta=0.1
    # expected 1/T_1 = 1 / (beta * pi * (eta^2 + 4)) at beta=100.0
    verify(
        TightBinding1D(; t=1.0, μ=0.0),
        NMRSpinRelaxationRate(),
        Infinite();
        route=:limiting_case,
        fetch_kw=(; beta=100.0, eta=0.1),
        independent=1.0 / (100.0 * π * (0.1^2 + 4.0)),
        agree_within=1e-5,
        refs=["Low-temperature regularized Korringa limit: 1/T_1 ~ T / (π * (η² + 4t²))"],
    )
end
# ── additional verification cards (#381 batch 7) ─────────────────────────
@testset "TightBinding1D — Energy/Infinite free fermion (#381 batch 7)" begin
    # Free spinless fermion at half-filling: e₀ = -(2t/π) sin(π·n_filling)
    # = -(2t/π) sin(π/2) = -2t/π (Ashcroft-Mermin 1976).
    for t in (0.5, 1.0, 2.0)
        verify(
            TightBinding1D(; t=t),
            Energy(:per_site),
            Infinite();
            route=:second_closed_form,
            independent=-2 * t / π,
            agree_within=1e-12,
            refs=["Ashcroft-Mermin 1976: 1D free spinless fermion half-filling e₀ = -2t/π"],
        )
    end
end

@testset "TightBinding1D — finite-T thermodynamics" begin
    # Finite-T integrals on the BZ.  All cards use verify(...) with an
    # independent analytic limit / textbook closed form; the implementation
    # is a black-box.  See test/util/verify.jl for the card schema.

    # high-T limit: ω(β → 0⁺; t, μ=0) → -log 2 / β.  Independent: textbook.
    for (t, β) in [(1.0, 1e-3), (2.5, 5e-4)]
        verify(
            QAtlas.TightBinding1D(; t=t, μ=0.0),
            QAtlas.FreeEnergy(),
            QAtlas.Infinite();
            route=:limiting_case,
            fetch_kw=(; beta=β),
            independent=(-log(2) / β),
            agree_within=abs(log(2) / β) * 2e-3,
            refs=[
                "Mahan, Many-Particle Physics §1.3: free-fermion β → 0⁺ limit ω → -T log 2 per site",
            ],
        )
    end

    # high-T limit: s(β → 0⁺) → log 2 per site (each mode half-occupied).
    for (t, μ, β) in [(1.0, 0.0, 1e-3), (1.0, 0.5, 1e-3), (2.0, -1.0, 5e-4)]
        verify(
            QAtlas.TightBinding1D(; t=t, μ=μ),
            QAtlas.ThermalEntropy(),
            QAtlas.Infinite();
            route=:limiting_case,
            fetch_kw=(; beta=β),
            independent=log(2),
            agree_within=log(2) * 4e-3,
            refs=[
                "Mahan, Many-Particle Physics §1.3: β → 0⁺ Fermi-Dirac entropy → log 2 per Bloch mode",
            ],
        )
    end

    # high-T limit: c_μ(β → 0⁺) → 0; bound by β² · (2t + |μ|)² (envelope of ε²).
    for (t, μ, β) in [(1.0, 0.0, 1e-2), (1.0, 0.5, 1e-2)]
        bound = (β * (2 * t + abs(μ)))^2
        verify(
            QAtlas.TightBinding1D(; t=t, μ=μ),
            QAtlas.SpecificHeat(),
            QAtlas.Infinite();
            route=:limiting_case,
            fetch_kw=(; beta=β),
            independent=0.0,
            agree_within=bound + 1e-12,
            refs=["Mahan, Many-Particle Physics §1.3: c_μ ~ β² · ⟨ε²⟩/4 → 0 as β → 0⁺"],
        )
    end

    # β → ∞ limit: ω(β=200) approaches T=0 grand-potential = -2/π at half-filling.
    # The T=0 value is independent (Ashcroft-Mermin Ch 9, free-fermion integral).
    verify(
        QAtlas.TightBinding1D(),
        QAtlas.FreeEnergy(),
        QAtlas.Infinite();
        route=:limiting_case,
        fetch_kw=(; beta=200.0),
        independent=-2 / π,
        agree_within=5e-3,
        refs=[
            "Ashcroft-Mermin (1976) Ch 9: half-filling 1D free-fermion E/N = -2/π = lim_{β→∞} ω(β)",
        ],
    )

    # Gibbs identity cross-check: s = β(u - ω).  Not a verify-card target
    # (relates three quantities, not one), kept as a plain @test sanity.
    @testset "Gibbs identity sanity (μ=0.3, β=2)" begin
        β = 2.0
        m = QAtlas.TightBinding1D(; t=1.0, μ=0.3)
        ω = QAtlas.fetch(m, QAtlas.FreeEnergy(), QAtlas.Infinite(); beta=β)
        s = QAtlas.fetch(m, QAtlas.ThermalEntropy(), QAtlas.Infinite(); beta=β)
        u = ω + s / β            # Gibbs: u = ω + s/β
        @test -2.3 ≤ u ≤ 2.3     # in-band internal energy
        cμ = QAtlas.fetch(m, QAtlas.SpecificHeat(), QAtlas.Infinite(); beta=β)
        @test cμ > 0             # gapless metal, strictly positive heat capacity
    end

    # DomainError on non-positive β — exception shape, kept as @test_throws.
    @testset "DomainError on β ≤ 0" begin
        @test_throws DomainError QAtlas.fetch(
            QAtlas.TightBinding1D(), QAtlas.FreeEnergy(), QAtlas.Infinite(); beta=0.0
        )
        @test_throws DomainError QAtlas.fetch(
            QAtlas.TightBinding1D(), QAtlas.ThermalEntropy(), QAtlas.Infinite(); beta=-1.0
        )
        @test_throws DomainError QAtlas.fetch(
            QAtlas.TightBinding1D(), QAtlas.SpecificHeat(), QAtlas.Infinite(); beta=0.0
        )
        @test_throws DomainError QAtlas.fetch(
            QAtlas.TightBinding1D(), QAtlas.NMRSpinRelaxationRate(), QAtlas.Infinite(); beta=0.0
        )
        @test_throws DomainError QAtlas.fetch(
            QAtlas.TightBinding1D(), QAtlas.NMRSpinRelaxationRate(), QAtlas.Infinite(); beta=1.0, eta=-0.1
        )
    end

    # ───────────────────────── NMRSpinRelaxationRate ──────────────────────────
    @testset "NMRSpinRelaxationRate — regularized 1/T_1" begin
        # High-T limit (beta -> 0): remains finite
        m = TightBinding1D(; t=1.0, μ=0.0)
        rate_high = QAtlas.fetch(m, NMRSpinRelaxationRate(), Infinite(); beta=1e-3, eta=0.1)
        @test rate_high > 0.0
        @test rate_high < 1.0

        # Low-T limit (beta -> infinity): matches regularized Korringa relation
        for eta_val in [0.05, 0.1, 0.2]
            expected_lim = 1.0 / (π * (eta_val^2 + 4.0))
            rate_low = QAtlas.fetch(m, NMRSpinRelaxationRate(), Infinite(); beta=100.0, eta=eta_val)
            @test isapprox(100.0 * rate_low, expected_lim; rtol=1e-3)
        end

        # Gapped insulator regime (|μ| > 2t): relaxation is exponentially suppressed at low-T
        m_gap = TightBinding1D(; t=1.0, μ=3.0)  # gap Δ = 1.0
        rate_gap_low = QAtlas.fetch(m_gap, NMRSpinRelaxationRate(), Infinite(); beta=10.0, eta=0.1)
        rate_gap_lower = QAtlas.fetch(m_gap, NMRSpinRelaxationRate(), Infinite(); beta=20.0, eta=0.1)
        @test rate_gap_lower < rate_gap_low * 0.1
    end
end
