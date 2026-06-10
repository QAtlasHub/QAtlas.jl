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
using QAtlas:
    TightBinding1D, Energy, MassGap, FermiVelocity, NMRSpinRelaxationRate, Infinite, fetch

# Independent reference for the η-regularized NMR rate: the same q-summed
# S(q,ω→0) golden-rule integral evaluated by a discrete midpoint k-mode sum over
# N modes — a DIFFERENT quadrature from the production nested QuadGK, converging
# to the continuum value as N→∞. Catches normalisation / Fermi-factor errors
# that a re-typed closed form cannot. (N=400 reproduces the QuadGK value to ~1e-9.)
_tb1d_nF(x) = x > 0 ? exp(-x) / (1 + exp(-x)) : 1 / (1 + exp(x))
function _tb1d_nmr_kmode_sum(t, μ, β, η; N=400)
    ks = [(n - 0.5) * π / N for n in 1:N]
    εs = [-2t * cos(k) - μ for k in ks]
    fs = _tb1d_nF.(β .* εs)
    s = 0.0
    for n in 1:N, m in 1:N
        s += fs[n] * (1 - fs[m]) * η / ((εs[n] - εs[m])^2 + η^2)
    end
    return s / (π * N^2)
end
# η-broadened particle–hole phase space (no Fermi factors); the high-T limit is
# 1/T₁(β→0) = ¼ · this, since f(1-f) → ¼ when every mode is half-filled.
function _tb1d_nmr_phasespace(t, μ, η; N=400)
    ks = [(n - 0.5) * π / N for n in 1:N]
    εs = [-2t * cos(k) - μ for k in ks]
    s = 0.0
    for n in 1:N, m in 1:N
        s += η / ((εs[n] - εs[m])^2 + η^2)
    end
    return s / (π * N^2)
end

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

    # NMR 1/T₁: the production nested-QuadGK integral vs an INDEPENDENT discrete
    # k-mode sum of the same q-summed S(q,ω→0) golden-rule integral (a different
    # quadrature) — this would catch any normalisation / Fermi-factor error.
    for (β, η) in ((1.0, 0.1), (0.5, 0.2), (2.0, 0.1))
        verify(
            TightBinding1D(; t=1.0, μ=0.0),
            NMRSpinRelaxationRate(),
            Infinite();
            route=:ed_finite_size,
            fetch_kw=(; beta=β, eta=η),
            independent=_tb1d_nmr_kmode_sum(1.0, 0.0, β, η; N=400),
            agree_within=1e-5,
            refs=[
                "Korringa 1950 (doi:10.1016/0031-8914(50)90105-4): free-fermion golden-rule 1/T₁ ∝ q-summed S(q,ω→0); cross-checked by a discrete N=400 k-mode sum → continuum.",
            ],
        )
    end
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
