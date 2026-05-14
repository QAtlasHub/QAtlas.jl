# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: TightBindingV1D (1D spinless fermion t-V chain).
#
# Targeted run (skips Pkg.test()):
#   julia --project=test test/models/quantum/misc/test_tight_binding_v1d.jl
#
# Phase 1 coverage (V = 0 free-fermion closed forms only):
#   • MassGap          = max(0, |μ| - 2t)
#   • FermiVelocity    = 2t · sin(arccos(-μ/(2t)))   for |μ| < 2t
#   • |μ| ≥ 2t         ⇒ FermiVelocity raises DomainError
#   • V ≠ 0            ⇒ MassGap / FermiVelocity raise DomainError (Phase 2)
#   • t ≤ 0            ⇒ constructor raises DomainError
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TightBindingV1D — V=0 free-fermion MassGap (Phase 1)" begin
    @test QAtlas.fetch(TightBindingV1D(), MassGap(), Infinite()) == 0.0      # default μ=0, gapless
    @test QAtlas.fetch(TightBindingV1D(; μ=1.5), MassGap(), Infinite()) == 0.0  # inside band
    @test QAtlas.fetch(TightBindingV1D(; μ=2.0), MassGap(), Infinite()) == 0.0  # band edge
    @test QAtlas.fetch(TightBindingV1D(; μ=3.0), MassGap(), Infinite()) == 1.0  # insulating
    @test QAtlas.fetch(TightBindingV1D(; μ=-5.0), MassGap(), Infinite()) == 3.0
    @test QAtlas.fetch(TightBindingV1D(; t=2.0, μ=5.0), MassGap(), Infinite()) == 1.0
end

@testset "TightBindingV1D — V=0 free-fermion FermiVelocity (Phase 1)" begin
    @test QAtlas.fetch(TightBindingV1D(), FermiVelocity(), Infinite()) ≈ 2.0          # μ=0 → v_F = 2t
    @test QAtlas.fetch(TightBindingV1D(; t=3.0), FermiVelocity(), Infinite()) ≈ 6.0   # v_F = 2·3
    @test QAtlas.fetch(TightBindingV1D(; μ=1.0), FermiVelocity(), Infinite()) ≈ sqrt(3)
    @test_throws DomainError QAtlas.fetch(
        TightBindingV1D(; μ=2.5), FermiVelocity(), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        TightBindingV1D(; μ=-3.0), FermiVelocity(), Infinite()
    )
end

@testset "TightBindingV1D — V≠0 throws DomainError (Phase 2 deferral)" begin
    @test_throws DomainError QAtlas.fetch(TightBindingV1D(; V=0.5), MassGap(), Infinite())
    @test_throws DomainError QAtlas.fetch(
        TightBindingV1D(; V=-1.0), FermiVelocity(), Infinite()
    )
end

@testset "TightBindingV1D — tiny V (1e-13) is non-zero (iszero strictness)" begin
    # Regression: V boundary must be exact iszero(V), not isapprox(...; atol=1e-12).
    # Any non-zero V puts us on the interacting JW-XXZ branch (Phase 2).
    @test_throws DomainError QAtlas.fetch(TightBindingV1D(; V=1e-13), MassGap(), Infinite())
    @test_throws DomainError QAtlas.fetch(
        TightBindingV1D(; V=-1e-13), FermiVelocity(), Infinite()
    )
end

@testset "TightBindingV1D — rejects t ≤ 0" begin
    @test_throws DomainError TightBindingV1D(; t=0.0)
    @test_throws DomainError TightBindingV1D(; t=-1.0)
end

@testset "TightBindingV1D — V=0 free-fermion Energy{:per_site} (Phase 1)" begin
    # μ=0 (half-filling, J=t=1): e₀ = -2/π
    e_hf = QAtlas.fetch(TightBindingV1D(), Energy(:per_site), Infinite())
    @test isapprox(e_hf, -2 / pi; atol=1e-12)
    @test isapprox(e_hf, -0.6366197723675814; atol=1e-12)

    # Band edges: empty (μ=-2t) → 0; filled (μ=+2t) → -μ = -2
    @test QAtlas.fetch(TightBindingV1D(; μ=-2.0), Energy(:per_site), Infinite()) == 0.0
    @test QAtlas.fetch(TightBindingV1D(; μ=2.0), Energy(:per_site), Infinite()) == -2.0
    # Deep insulating limits
    @test QAtlas.fetch(TightBindingV1D(; μ=-5.0), Energy(:per_site), Infinite()) == 0.0
    @test QAtlas.fetch(TightBindingV1D(; μ=5.0), Energy(:per_site), Infinite()) == -5.0

    # t-linearity at μ=0: e₀(t, 0) = -2t/π
    @test isapprox(
        QAtlas.fetch(TightBindingV1D(; t=3.0), Energy(:per_site), Infinite()),
        -2 * 3 / pi;
        atol=1e-12,
    )

    # Continuity across band edges (one-sided)
    e_just_inside = QAtlas.fetch(
        TightBindingV1D(; μ=2.0 - 1e-9), Energy(:per_site), Infinite()
    )
    @test isapprox(e_just_inside, -2.0; atol=1e-6)
    e_just_empty = QAtlas.fetch(
        TightBindingV1D(; μ=-2.0 + 1e-9), Energy(:per_site), Infinite()
    )
    @test isapprox(e_just_empty, 0.0; atol=1e-6)
end

@testset "TightBindingV1D — Energy{:per_site} V≠0 / tiny-V DomainError (Phase 2 gate)" begin
    @test_throws DomainError QAtlas.fetch(
        TightBindingV1D(; V=0.5), Energy(:per_site), Infinite()
    )
    # iszero(V) strictness regression
    @test_throws DomainError QAtlas.fetch(
        TightBindingV1D(; V=1e-13), Energy(:per_site), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        TightBindingV1D(; V=-1e-13), Energy(:per_site), Infinite()
    )
end
