using QAtlas, Test

@testset "TFIM — CriticalExponents = 2D Ising Onsager (Phase 2)" begin
    m = TFIM()
    exp = QAtlas.fetch(m, CriticalExponents(), Infinite())
    @test exp.β == 1 // 8
    @test exp.γ == 7 // 4
    @test exp.δ == 15
    @test exp.ν == 1
    @test exp.α == 0
    @test exp.η == 1 // 4
    # Delegation invariant
    @test exp == QAtlas.fetch(QAtlas.Universality(:Ising), CriticalExponents(); d=2)
    # Rushbrooke α + 2β + γ = 2
    @test exp.α + 2 * exp.β + exp.γ == 2
    # Widom γ = β(δ − 1)
    @test exp.γ == exp.β * (exp.δ - 1)
    # Fisher η = 2 − γ/ν (TFIM: 2 − 7/4 = 1/4)
    @test exp.η == 2 - exp.γ // exp.ν
end
