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
