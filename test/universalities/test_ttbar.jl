using QAtlas, Test

@testset "TTbar — UV CentralCharge invariance (Phase 1)" begin
    # Default: c = 1, λ = 0 → undeformed CFT, c = 1
    @test QAtlas.fetch(TTbar(), CentralCharge(), Infinite()) == 1.0
    # c is preserved at any λ ≠ 0 (TT̄ is irrelevant)
    for λ in (-1.0, 0.0, 0.5, 2.0, -10.0)
        @test QAtlas.fetch(TTbar(; c=1.0, λ=λ), CentralCharge(), Infinite()) == 1.0
    end
    # Various seed UV central charges
    for c in (0.5, 1.0, 2.0, 14.0)
        @test QAtlas.fetch(TTbar(; c=c, λ=0.7), CentralCharge(), Infinite()) == c
    end
end

@testset "TTbar — rejects c ≤ 0 (Phase 1)" begin
    @test_throws DomainError TTbar(; c=0.0)
    @test_throws DomainError TTbar(; c=-1.0)
    m = TTbar(; c=1.0)
    @test_throws DomainError QAtlas.fetch(m, CentralCharge(), Infinite(); c=0.0)
end
