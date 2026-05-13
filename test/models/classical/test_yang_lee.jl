# ─────────────────────────────────────────────────────────────────────────────
# Test: YangLee — c = -22/5 and Kac primaries via MinimalModel(5, 2).
#
# Verifies (Phase 1):
#   * Central charge exactly -22//5 (Rational, machine-precision agreement).
#   * Delegation invariant: equals MinimalModel(5, 2)'s c.
#   * Kac primaries h_{1,1} = 0 and h_{1,2} = -1/5 (the famous Yang-Lee
#     negative-dimension primary), and their Kac-symmetric duals
#     h_{1,4} = 0, h_{1,3} = -1/5.
#   * Index range guards: r ∈ [1, 1], s ∈ [1, 4] raise DomainError.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "YangLee — CentralCharge c = -22/5 (Phase 1, M(5,2))" begin
    c = QAtlas.fetch(YangLee(), CentralCharge(), Infinite())
    @test c == -22 // 5
    @test c isa Rational
    @test c isa Rational{Int}
    # Delegation invariant
    @test c == QAtlas.fetch(QAtlas.MinimalModel(5, 2), CentralCharge())
end

@testset "YangLee — ConformalWeights (Phase 1)" begin
    m = YangLee()
    # Identity
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); r=1, s=1) == 0
    # The famous negative-dimension Yang-Lee primary
    h12 = QAtlas.fetch(m, ConformalWeights(), Infinite(); r=1, s=2)
    @test h12 == -1 // 5
    @test h12 isa Rational
    # Kac symmetry (1, s) ↔ (1, p - s) for p = 5
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); r=1, s=3) == -1 // 5
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); r=1, s=4) == 0
    # Delegation invariant on each Kac primary (not just c)
    for s in 1:4
        @test QAtlas.fetch(m, ConformalWeights(), Infinite(); r=1, s=s) ==
            QAtlas.fetch(QAtlas.MinimalModel(5, 2), ConformalWeights(); r=1, s=s)
    end
end

@testset "YangLee — index range guards (Phase 1)" begin
    m = YangLee()
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); r=0, s=1)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); r=2, s=1)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); r=1, s=0)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); r=1, s=5)
end
