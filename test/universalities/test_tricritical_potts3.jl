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
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); r=1, s=1) == 0
    # Energy operator ε: h_{1,2} = 1/7
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); r=1, s=2) == 1 // 7
    # h_{2,1} = 3/8
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); r=2, s=1) == 3 // 8
    # h_{2,2} = 1/56
    @test QAtlas.fetch(m, ConformalWeights(), Infinite(); r=2, s=2) == 1 // 56
    # Delegation invariant: matches MinimalModel(7, 6) exactly
    for (r, s) in [(1, 1), (1, 2), (2, 1), (2, 2), (3, 3), (5, 6)]
        @test QAtlas.fetch(m, ConformalWeights(), Infinite(); r=r, s=s) ==
            QAtlas.fetch(QAtlas.MinimalModel(7, 6), ConformalWeights(); r=r, s=s)
    end
end

@testset "TricriticalPotts3 — ConformalWeights index range (Phase 2)" begin
    m = TricriticalPotts3()
    # r ∈ [1, 5], s ∈ [1, 6] for M(7, 6)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); r=0, s=1)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); r=6, s=1)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); r=1, s=0)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); r=1, s=7)
end

@testset "TricriticalPotts3 — PrimaryFields delegation (Phase 2)" begin
    m = TricriticalPotts3()
    pf_tp = QAtlas.fetch(m, PrimaryFields(), Infinite())
    pf_mm = QAtlas.fetch(QAtlas.MinimalModel(7, 6), PrimaryFields())
    @test pf_tp == pf_mm
    # M(7, 6) has (p_prime - 1)*(p - 1)/2 = 5*6/2 = 15 independent primaries
    @test length(pf_tp) == 15
    # Identity primary appears with h = 0
    @test any(x -> x.h == 0, pf_tp)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "TricriticalPotts3 — verification cards" begin
    # Tricritical 3-state Potts = minimal model M(7,6): c = 6/7
    # (independent Kac/minimal-model formula 1 - 6(p-q)²/(pq)).
    verify(
        TricriticalPotts3(),
        CentralCharge(),
        Infinite();
        route=:second_closed_form,
        independent=1 - 6 * (7 - 6)^2 / (7 * 6),
        agree_within=1e-12,
        refs=["M(7,6): c = 1 - 6(p-q)²/(pq) = 6/7"],
    )
end
# ── additional verification cards (#381 batch) ─────────────────────────────
@testset "TricriticalPotts3 — additional verification cards (#381 batch)" begin
    # Identity primary h_{1,1} = 0 (Kac formula, M(7,6)).
    verify(
        TricriticalPotts3(),
        ConformalWeights(),
        Infinite();
        route=:second_closed_form,
        independent=0.0,
        agree_within=1e-14,
        refs=["BPZ 1984; Andrews-Baxter-Forrester 1984: identity h_{1,1} = 0 in M(7,6)"],
        fetch_kw=(; r=1, s=1),
    )
end
