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

    # Explicit central charge of the universality class — c = 1/2 (Onsager).
    @test QAtlas.fetch(QAtlas.Universality(:Ising), CentralCharge(); d=2) == 1 // 2

    # Hyperscaling relations
    @test exp.α + 2 * exp.β + exp.γ == 2          # Rushbrooke
    @test exp.γ == exp.β * (exp.δ - 1)            # Widom
    @test exp.η == 2 - exp.γ // exp.ν             # Fisher
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "IsingSquare critical — verification cards" begin
    # CriticalExponents() returns a NamedTuple (non-scalar; covered by
    # the delegation tests above).  The scalar Onsager Tc anchors the
    # same critical point with an independent closed form.
    verify(
        IsingSquare(; J=1.0),
        CriticalTemperature(),
        Infinite();
        route=:second_closed_form,
        independent=2.0 / log(1 + sqrt(2)),
        agree_within=1e-10,
        refs=["Onsager 1944: Tc = 2J / log(1+√2) ≈ 2.269185 (β=1/8 universality anchor)"],
    )
end
