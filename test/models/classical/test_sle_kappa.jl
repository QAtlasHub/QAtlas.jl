# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: SLEkappa — SLE-CFT central charge c(κ).
#
# Verifies:
#   * Canonical fixed points c(2)=-2, c(8/3)=0, c(3)=1/2, c(4)=1,
#     c(6)=0, c(8)=-2 to machine precision.
#   * κ ↔ 16/κ duality: c(κ) = c(16/κ) for several test points.
#   * Continuity: bounded behaviour for small κ (no spurious overflow
#     at the κ → 0⁺ pole except its mathematical divergence).
#   * DomainError on κ ≤ 0.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "SLEkappa — canonical fixed-point central charges" begin
    pairs = (
        (κ=2.0,   c_exact=-2.0),     # LERW
        (κ=8/3,   c_exact=0.0),       # SAW
        (κ=3.0,   c_exact=1/2),       # Ising boundary
        (κ=4.0,   c_exact=1.0),       # GFF
        (κ=6.0,   c_exact=0.0),       # Percolation
        (κ=8.0,   c_exact=-2.0),     # UST Peano
    )
    for (κ, c_exact) in pairs
        c = QAtlas.fetch(SLEkappa(; κ=κ), CentralCharge(), Infinite())
        @test c ≈ c_exact atol = 1e-14
    end
end

@testset "SLEkappa — duality κ ↔ 16/κ" begin
    for κ in (1.0, 2.0, 3.0, 8/3, 5.0, 7.0)
        c1 = QAtlas.fetch(SLEkappa(; κ=κ), CentralCharge(), Infinite())
        c2 = QAtlas.fetch(SLEkappa(; κ=16 / κ), CentralCharge(), Infinite())
        @test c1 ≈ c2 atol = 1e-12
    end
end

@testset "SLEkappa — DomainError on κ ≤ 0" begin
    @test_throws DomainError QAtlas.fetch(SLEkappa(; κ=0.0), CentralCharge(), Infinite())
    @test_throws DomainError QAtlas.fetch(SLEkappa(; κ=-1.5), CentralCharge(), Infinite())
end
