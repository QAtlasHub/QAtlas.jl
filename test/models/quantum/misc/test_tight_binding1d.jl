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
#     t-linearity, DomainError in the gapped phase.
#   • Constructor — rejects t ≤ 0 with DomainError.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test
using QAtlas: TightBinding1D, Energy, MassGap, FermiVelocity, Infinite, fetch

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

        # No Fermi surface in the insulating regime — both signs of μ.
        @test_throws DomainError QAtlas.fetch(
            TightBinding1D(; t=1.0, μ=2.5), FermiVelocity(), Infinite()
        )
        @test_throws DomainError QAtlas.fetch(
            TightBinding1D(; t=1.0, μ=-3.0), FermiVelocity(), Infinite()
        )
    end

    # ────────────────────────── Constructor validation ─────────────────────────
    @testset "Constructor rejects t ≤ 0" begin
        @test_throws DomainError TightBinding1D(; t=0.0)
        @test_throws DomainError TightBinding1D(; t=-1.0)
    end
end
