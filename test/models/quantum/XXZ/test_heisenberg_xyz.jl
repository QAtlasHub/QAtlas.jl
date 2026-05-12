# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: HeisenbergXYZ — axis-aligned XXZ-delegation.
#
# Verifies:
#   * Jx = Jy = J, Jz = J (isotropic AF Heisenberg):
#         E0/N = J (1/4 - log 2)   (Hulthén 1938 via XXZ1D Δ = 1 path)
#   * Jx = Jy = J, Jz = 0 (XX free fermion): E0/N = -J/π
#   * Jx = Jy = J, Jz = -J: E0/N = -J/4 (FM saturation)
#   * Jx = Jy = J, Jz = J/2: matches XXZ1D Yang-Yang Δ = 1/2 → -3J/8
#   * Jx ≠ Jy: DomainError (general XYZ deferred to Phase 2)
#   * Jx = Jy = 0: DomainError (Ising-like reduction)
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "HeisenbergXYZ — isotropic Heisenberg AF limit (Hulthén)" begin
    m = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=1.0)
    E0 = QAtlas.fetch(m, Energy(:per_site), Infinite())
    @test E0 ≈ 0.25 - log(2.0) atol = 1e-12
end

@testset "HeisenbergXYZ — XX free-fermion limit" begin
    m = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=0.0)
    E0 = QAtlas.fetch(m, Energy(:per_site), Infinite())
    @test E0 ≈ -1 / π atol = 1e-12
end

@testset "HeisenbergXYZ — isotropic FM saturated limit" begin
    m = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=-1.0)
    E0 = QAtlas.fetch(m, Energy(:per_site), Infinite())
    @test E0 ≈ -1 / 4 atol = 1e-12
end

@testset "HeisenbergXYZ — Yang-Yang single integral at Δ = 1/2" begin
    m = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=0.5)
    E0 = QAtlas.fetch(m, Energy(:per_site), Infinite())
    @test E0 ≈ -3 / 8 atol = 1e-10
end

@testset "HeisenbergXYZ — J-scaling delegated to XXZ1D" begin
    m1 = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=0.5)
    m2 = HeisenbergXYZ(; Jx=2.5, Jy=2.5, Jz=1.25)
    e1 = QAtlas.fetch(m1, Energy(:per_site), Infinite())
    e2 = QAtlas.fetch(m2, Energy(:per_site), Infinite())
    # XXZ1D returns J × (Δ-only function), so e2 = 2.5 e1.
    @test e2 ≈ 2.5 * e1 atol = 1e-10
end

@testset "HeisenbergXYZ — Jx ≠ Jy raises DomainError" begin
    @test_throws DomainError QAtlas.fetch(
        HeisenbergXYZ(; Jx=1.0, Jy=0.5, Jz=0.3), Energy(:per_site), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        HeisenbergXYZ(; Jx=1.0, Jy=2.0, Jz=1.0), Energy(:per_site), Infinite()
    )
end

@testset "HeisenbergXYZ — Jx = Jy = 0 raises DomainError" begin
    @test_throws DomainError QAtlas.fetch(
        HeisenbergXYZ(; Jx=0.0, Jy=0.0, Jz=1.0), Energy(:per_site), Infinite()
    )
end
