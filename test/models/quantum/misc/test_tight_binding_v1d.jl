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

@testset "TightBindingV1D — rejects t ≤ 0" begin
    @test_throws DomainError TightBindingV1D(; t=0.0)
    @test_throws DomainError TightBindingV1D(; t=-1.0)
end
