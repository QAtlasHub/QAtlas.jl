# test/models/classical/test_random_bond_ising_2d.jl
#
# RandomBondIsing2D — Phase 1 tests (p = 1 pure FM → 2D Ising c = 1/2).

using QAtlas, Test

@testset "RandomBondIsing2D — p=1 FM critical c=1/2 (Phase 1)" begin
    c = QAtlas.fetch(RandomBondIsing2D(; J=1.0, p=1.0), CentralCharge(), Infinite())
    @test c == 1 // 2
    # Default constructor IS the pure FM point.
    @test QAtlas.fetch(RandomBondIsing2D(), CentralCharge(), Infinite()) == 1 // 2
    # Delegation invariant.
    @test c == QAtlas.fetch(QAtlas.MinimalModel(4, 3), CentralCharge())
end

@testset "RandomBondIsing2D — Nishimori/disordered p ≠ 1 throws DomainError (Phase 2)" begin
    @test_throws DomainError QAtlas.fetch(
        RandomBondIsing2D(; p=0.5), CentralCharge(), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        RandomBondIsing2D(; p=0.1093), CentralCharge(), Infinite()
    )  # Nishimori MC point
    @test_throws DomainError QAtlas.fetch(
        RandomBondIsing2D(; p=0.0), CentralCharge(), Infinite()
    )  # all anti-FM
end

@testset "RandomBondIsing2D — rejects J ≤ 0 or p ∉ [0,1] (Phase 1)" begin
    @test_throws DomainError RandomBondIsing2D(; J=0.0)
    @test_throws DomainError RandomBondIsing2D(; J=-1.0)
    @test_throws DomainError RandomBondIsing2D(; p=-0.1)
    @test_throws DomainError RandomBondIsing2D(; p=1.5)
end
