# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: LiouvilleCFT — c = 1 + 6(b + 1/b)².
#
# Verifies:
#   * Self-dual point b = 1: c = 25 exactly.
#   * b ↔ 1/b duality: c(b) = c(1/b) for several test points.
#   * c → ∞ as b → 0⁺ and b → ∞ (large-coupling divergence both ways).
#   * DomainError on b ≤ 0.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "LiouvilleCFT — self-dual point b = 1 gives c = 25" begin
    @test QAtlas.fetch(LiouvilleCFT(; b=1.0), CentralCharge(), Infinite()) ≈ 25.0 atol = 1e-14
end

@testset "LiouvilleCFT — b ↔ 1/b duality" begin
    for b in (0.5, 0.7, 1.3, 2.0, 5.0)
        c1 = QAtlas.fetch(LiouvilleCFT(; b=b), CentralCharge(), Infinite())
        c2 = QAtlas.fetch(LiouvilleCFT(; b=1 / b), CentralCharge(), Infinite())
        @test c1 ≈ c2 atol = 1e-12
    end
end

@testset "LiouvilleCFT — small-b and large-b growth" begin
    c_small = QAtlas.fetch(LiouvilleCFT(; b=0.01), CentralCharge(), Infinite())
    c_large = QAtlas.fetch(LiouvilleCFT(; b=100.0), CentralCharge(), Infinite())
    @test c_small > 1e3
    @test c_large > 1e3
    @test c_small ≈ c_large atol = 1e-8     # exact duality at b = 1/100 ↔ b = 100
end

@testset "LiouvilleCFT — DomainError on b ≤ 0" begin
    @test_throws DomainError QAtlas.fetch(
        LiouvilleCFT(; b=0.0), CentralCharge(), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        LiouvilleCFT(; b=-1.0), CentralCharge(), Infinite()
    )
end
