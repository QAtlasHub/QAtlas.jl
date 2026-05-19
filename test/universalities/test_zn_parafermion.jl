# ─────────────────────────────────────────────────────────────────────────────
# Universalities test: ZnParafermion — c = 2(n-1)/(n+2) (Phase 1).
#
# Verifies:
#   * Special values at n ∈ {2, 3, 4, 5, 6} (exact Rational).
#   * Default constructor selects Z_3 parafermion (c = 4/5).
#   * n = 2 matches MinimalModel(4, 3) (Ising delegation invariant).
#   * Asymptote c(n→∞) → 2 (large-n sanity check).
#   * Domain error for n < 2.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "ZnParafermion — CentralCharge c = 2(n-1)/(n+2) (Phase 1)" begin
    # Special values
    @test QAtlas.fetch(ZnParafermion(; n=2), CentralCharge(), Infinite()) == 1 // 2   # Ising
    @test QAtlas.fetch(ZnParafermion(; n=3), CentralCharge(), Infinite()) == 4 // 5   # 3-state Potts (default)
    @test QAtlas.fetch(ZnParafermion(), CentralCharge(), Infinite()) == 4 // 5         # default = Z_3
    @test QAtlas.fetch(ZnParafermion(; n=4), CentralCharge(), Infinite()) == 1 // 1   # = 1 free boson
    @test QAtlas.fetch(ZnParafermion(; n=5), CentralCharge(), Infinite()) == 8 // 7
    @test QAtlas.fetch(ZnParafermion(; n=6), CentralCharge(), Infinite()) == 10 // 8
    # n → ∞ approaches c = 2
    @test QAtlas.fetch(ZnParafermion(; n=1000), CentralCharge(), Infinite()) ≈ 2.0 atol=0.01
end

@testset "ZnParafermion — Ising at n=2 matches MinimalModel(4,3) (Phase 1)" begin
    @test QAtlas.fetch(ZnParafermion(; n=2), CentralCharge(), Infinite()) ==
        QAtlas.fetch(QAtlas.MinimalModel(4, 3), CentralCharge())
end

@testset "ZnParafermion — Z_3 Potts at n=3 matches MinimalModel(6,5) (Phase 1)" begin
    # 3-state Potts critical CFT: Z_3 parafermion = M(6, 5), both c = 4/5.
    @test QAtlas.fetch(ZnParafermion(; n=3), CentralCharge(), Infinite()) ==
        QAtlas.fetch(QAtlas.MinimalModel(6, 5), CentralCharge())
end

@testset "ZnParafermion — rejects n < 2 (Phase 1)" begin
    @test_throws DomainError ZnParafermion(; n=1)
    @test_throws DomainError ZnParafermion(; n=0)
    @test_throws DomainError ZnParafermion(; n=-3)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "ZnParafermion — verification cards" begin
    # Fateev-Zamolodchikov Z_n parafermion CFT: c = 2(n-1)/(n+2)
    # (independent closed form).
    for n in (2, 3, 4, 5)
        verify(
            ZnParafermion(; n=n),
            CentralCharge(),
            Infinite();
            route=:second_closed_form,
            independent=2 * (n - 1) // (n + 2),
            agree_within=1e-12,
            refs=["Fateev-Zamolodchikov: Z_n parafermion c = 2(n-1)/(n+2)"],
        )
    end
end
