# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: TricriticalPotts3 — c = 6/7 via MinimalModel(6,7).
#
# Verifies:
#   * Central charge exactly 6//7 (Rational, machine-precision agreement).
#   * Result equals MinimalModel(7, 6)'s c (delegation invariant).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TricriticalPotts3 — c = 6/7 exact" begin
    c = QAtlas.fetch(TricriticalPotts3(), CentralCharge(), Infinite())
    @test c == 6 // 7
end

@testset "TricriticalPotts3 — equals MinimalModel(7, 6)" begin
    c_tp = QAtlas.fetch(TricriticalPotts3(), CentralCharge(), Infinite())
    c_mm = QAtlas.fetch(QAtlas.MinimalModel(7, 6), CentralCharge())
    @test c_tp == c_mm
end

@testset "TricriticalPotts3 — ConformalWeights delegation (Phase 2)" begin
    m = TricriticalPotts3()
    # Identity has h = 0
    @test QAtlas.fetch(m, ConformalWeights(); r=1, s=1) == 0
    # Energy operator ε: h_{1,2} = 1/7
    @test QAtlas.fetch(m, ConformalWeights(); r=1, s=2) == 1 // 7
    # h_{2,1} = 3/8
    @test QAtlas.fetch(m, ConformalWeights(); r=2, s=1) == 3 // 8
    # h_{2,2} = 1/56
    @test QAtlas.fetch(m, ConformalWeights(); r=2, s=2) == 1 // 56
    # Delegation invariant: matches MinimalModel(7, 6) exactly
    for (r, s) in [(1, 1), (1, 2), (2, 1), (2, 2), (3, 3), (5, 6)]
        @test QAtlas.fetch(m, ConformalWeights(); r=r, s=s) ==
              QAtlas.fetch(QAtlas.MinimalModel(7, 6), ConformalWeights(); r=r, s=s)
    end
end

@testset "TricriticalPotts3 — ConformalWeights index range (Phase 2)" begin
    m = TricriticalPotts3()
    # r ∈ [1, 5], s ∈ [1, 6] for M(7, 6)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(); r=0, s=1)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(); r=6, s=1)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(); r=1, s=0)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(); r=1, s=7)
end

@testset "TricriticalPotts3 — PrimaryFields delegation (Phase 2)" begin
    m = TricriticalPotts3()
    pf_tp = QAtlas.fetch(m, PrimaryFields())
    pf_mm = QAtlas.fetch(QAtlas.MinimalModel(7, 6), PrimaryFields())
    @test pf_tp == pf_mm
    # M(7, 6) has (p_prime - 1)*(p - 1)/2 = 5*6/2 = 15 independent primaries
    @test length(pf_tp) == 15
    # Identity primary appears with h = 0
    @test any(x -> x.h == 0, pf_tp)
end
