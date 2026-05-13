# ─────────────────────────────────────────────────────────────────────────────
# Test: TricriticalIsing — M(5,4) unitary minimal model, c = 7/10.
#
# Verifies (Phase 1):
#   * c = 7/10 exact (Rational), matches MinimalModel(5, 4) delegation.
#   * Kac weights h_{r,s} for the six famous primaries (σ, σ', ε, ε', ε'', 1).
#   * PrimaryFields enumeration has 6 entries and matches MinimalModel.
#   * Out-of-range r, s raise DomainError (forwarded from MinimalModel).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TricriticalIsing — CentralCharge c = 7/10 (Phase 1, M(5,4))" begin
    c = QAtlas.fetch(TricriticalIsing(), CentralCharge(), Infinite())
    @test c == 7 // 10
    # Delegation invariant
    @test c == QAtlas.fetch(QAtlas.MinimalModel(5, 4), CentralCharge())
end

@testset "TricriticalIsing — ConformalWeights (Phase 1)" begin
    m = TricriticalIsing()
    # MinimalModel API uses (r, s) with 1 ≤ r ≤ p_prime−1 = 3,
    # 1 ≤ s ≤ p−1 = 4, h = ((p r − p_prime s)² − (p − p_prime)²) / (4 p p_prime)
    # = ((5 r − 4 s)² − 1) / 80.  Famous primaries below are quoted
    # in the literature (r̃, s̃) Kac convention with the API's
    # (r, s) = (s̃, r̃) substitution.
    # Identity — Kac (1,1)
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); r=1, s=1) == 0
    # σ (spin) — Kac (2,2) (symmetric in r ↔ s under p ↔ p_prime swap)
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); r=2, s=2) == 3 // 80
    # σ′ (subleading spin) — Kac literature (1,2) → API (r=2, s=1)
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); r=2, s=1) == 7 // 16
    # ε (energy) — Kac literature (2,1) → API (r=1, s=2)
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); r=1, s=2) == 1 // 10
    # ε′ (vacancy) — Kac literature (2,3) → API (r=3, s=2)
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); r=3, s=2) == 3 // 5
    # ε″ (irrelevant) — Kac literature (1,3) → API (r=3, s=1)
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); r=3, s=1) == 3 // 2
end

@testset "TricriticalIsing — PrimaryFields (Phase 1)" begin
    m = TricriticalIsing()
    pf = QAtlas.fetch(m, PrimaryFields(), Infinite())
    @test length(pf) == 6   # (p−1)(p_prime−1)/2 = 4·3/2 = 6
    @test pf == QAtlas.fetch(QAtlas.MinimalModel(5, 4), PrimaryFields())
end

@testset "TricriticalIsing — ConformalWeights index range" begin
    m = TricriticalIsing()
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); r=0, s=1)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); r=4, s=1)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); r=1, s=0)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); r=1, s=5)
end
