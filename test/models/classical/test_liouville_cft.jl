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
    @test QAtlas.fetch(LiouvilleCFT(; b=1.0), CentralCharge(), Infinite()) ≈ 25.0 atol =
        1e-14
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

# ─────────────────────────────────────────────────────────────────────────────
# Phase 2: ConformalWeights — Δ_α = α(Q − α).
# ─────────────────────────────────────────────────────────────────────────────

@testset "LiouvilleCFT — ConformalWeights identity (α = 0 ⇒ Δ = 0)" begin
    for b in (0.5, 1.0, 2.0, 3.7)
        @test QAtlas.fetch(LiouvilleCFT(; b=b), ConformalWeights(), Infinite(); α=0.0) ≈ 0.0 atol =
            1e-14
    end
end

@testset "LiouvilleCFT — ConformalWeights at b = 1 (Q = 2)" begin
    m = LiouvilleCFT(; b=1.0)
    # α = 1 = b = Q/2 (triply degenerate at the self-dual point)
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); α=1.0) ≈ 1.0 atol = 1e-14
    # α = Q = 2 → Δ = 0 (reflection of identity)
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); α=2.0) ≈ 0.0 atol = 1e-14
end

@testset "LiouvilleCFT — ConformalWeights at b = 2 (Q = 2.5, non-degenerate)" begin
    m = LiouvilleCFT(; b=2.0)
    # α = b = 2 → Δ = 2·(2.5 − 2) = 1 (degenerate screening V_b)
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); α=2.0) ≈ 1.0 atol = 1e-14
    # α = 1/b = 0.5 → Δ = 0.5·(2.5 − 0.5) = 1 (dual screening V_{1/b})
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); α=0.5) ≈ 1.0 atol = 1e-14
    # α = Q/2 = 1.25 → Δ = Q²/4 = 1.5625 (Seiberg bound)
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); α=1.25) ≈ 1.5625 atol = 1e-14
end

@testset "LiouvilleCFT — ConformalWeights reflection symmetry Δ_α = Δ_{Q−α}" begin
    b = 0.7
    Q = b + 1 / b
    m = LiouvilleCFT(; b=b)
    for α in (-0.3, 0.0, 0.25, 0.5 * Q, 0.9, 1.4, Q, Q + 0.5)
        Δ_α = QAtlas.fetch(m, ConformalWeights(), Infinite(); α=α)
        Δ_refl = QAtlas.fetch(m, ConformalWeights(), Infinite(); α=Q - α)
        @test Δ_α ≈ Δ_refl atol = 1e-12
    end
end

@testset "LiouvilleCFT — ConformalWeights DomainError on b ≤ 0" begin
    @test_throws DomainError QAtlas.fetch(
        LiouvilleCFT(; b=0.0), ConformalWeights(), Infinite(); α=0.5
    )
    @test_throws DomainError QAtlas.fetch(
        LiouvilleCFT(; b=-1.0), ConformalWeights(), Infinite(); α=0.5
    )
end
