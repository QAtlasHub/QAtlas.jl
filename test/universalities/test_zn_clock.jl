using QAtlas, Test

@testset "ZnClock — n=2 Ising universality (Phase 1)" begin
    @test QAtlas.fetch(ZnClock(; n=2), CentralCharge(), Infinite()) == 1 // 2
    # Default constructor
    @test QAtlas.fetch(ZnClock(), CentralCharge(), Infinite()) == 1 // 2
    # Delegation invariant
    @test QAtlas.fetch(ZnClock(; n=2), CentralCharge(), Infinite()) ==
        QAtlas.fetch(QAtlas.MinimalModel(4, 3), CentralCharge())
end

@testset "ZnClock — n=3 Potts universality (Phase 1)" begin
    @test QAtlas.fetch(ZnClock(; n=3), CentralCharge(), Infinite()) == 4 // 5
    @test QAtlas.fetch(ZnClock(; n=3), CentralCharge(), Infinite()) ==
        QAtlas.fetch(QAtlas.MinimalModel(6, 5), CentralCharge())
end

@testset "ZnClock — n ≥ 4 throws DomainError (Phase 2 deferral)" begin
    @test_throws DomainError QAtlas.fetch(ZnClock(; n=4), CentralCharge(), Infinite())
    @test_throws DomainError QAtlas.fetch(ZnClock(; n=5), CentralCharge(), Infinite())
    @test_throws DomainError QAtlas.fetch(ZnClock(; n=6), CentralCharge(), Infinite())
    @test_throws DomainError QAtlas.fetch(ZnClock(; n=100), CentralCharge(), Infinite())
end

@testset "ZnClock — rejects n < 2 (Phase 1)" begin
    @test_throws DomainError ZnClock(; n=1)
    @test_throws DomainError ZnClock(; n=0)
    @test_throws DomainError ZnClock(; n=-1)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "ZnClock — verification cards" begin
    # n=2 clock = Ising = M(4,3): c = 1 - 6(p-q)²/(pq) = 1/2
    verify(
        ZnClock(; n=2),
        CentralCharge(),
        Infinite();
        route=:second_closed_form,
        independent=1 - 6 * (4 - 3)^2 // (4 * 3),
        agree_within=1e-12,
        refs=["n=2 clock = Ising = M(4,3): c = 1/2"],
    )
    # n=3 clock = 3-state Potts = M(6,5): c = 4/5
    verify(
        ZnClock(; n=3),
        CentralCharge(),
        Infinite();
        route=:second_closed_form,
        independent=1 - 6 * (6 - 5)^2 // (6 * 5),
        agree_within=1e-12,
        refs=["n=3 clock = 3-state Potts = M(6,5): c = 4/5"],
    )
end
