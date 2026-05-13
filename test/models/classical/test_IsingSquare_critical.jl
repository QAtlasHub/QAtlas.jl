# Phase 2: IsingSquare critical exponents — delegation to Universality(:Ising) d=2.
# Verifies Onsager 1944 exact 2D Ising universality and the hyperscaling
# relations (Rushbrooke, Widom, Fisher).

using QAtlas, Test

@testset "IsingSquare — CriticalExponents = Onsager 2D Ising (Phase 2)" begin
    m = IsingSquare()
    exp = QAtlas.fetch(m, CriticalExponents(), Infinite())

    # Onsager 1944 exact exponents
    @test exp.α == 0
    @test exp.β == 1 // 8
    @test exp.γ == 7 // 4
    @test exp.δ == 15
    @test exp.ν == 1
    @test exp.η == 1 // 4

    # Delegation invariant: identical payload as the universality entry.
    @test exp == QAtlas.fetch(QAtlas.Universality(:Ising), CriticalExponents(); d=2)

    # Hyperscaling relations
    @test exp.α + 2 * exp.β + exp.γ == 2          # Rushbrooke
    @test exp.γ == exp.β * (exp.δ - 1)            # Widom
    @test exp.η == 2 - exp.γ // exp.ν             # Fisher
end
