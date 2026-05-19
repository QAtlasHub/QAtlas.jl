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

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "YangLee — verification cards" begin
    # M(5,2) minimal model: c = 1 - 6(p-p')²/(p p') with p=5, p'=2
    verify(
        YangLee(),
        CentralCharge(),
        Infinite();
        route=:second_closed_form,
        independent=1 - 6 * (5 - 2)^2 / (5 * 2),
        agree_within=1e-12,
        refs=["Yang-Lee edge = M(5,2): c = 1 - 6(p-p')²/(pp') = -22/5"],
    )

    # Yang-Lee edge primary conformal weight h = -1/5 (Cardy 1985)
    verify(
        YangLee(),
        ConformalWeights(),
        Infinite();
        route=:literature_value,
        fetch_kw=(; r=1, s=2),
        independent=-1 / 5,
        agree_within=1e-12,
        refs=["Cardy 1985: Yang-Lee edge singularity primary h_{1,2} = -1/5"],
    )
end
