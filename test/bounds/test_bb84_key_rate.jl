# Shor-Preskill BB84 secret-key rate — a *lower* bound (Bound{:QuantumInformation}),
# resolved from the deferred dump.  R(e) = 1 - 2 H₂(e), achievable key fraction.
# (verify_bound is in scope via the suite's global include of test/util/verify.jl.)

using QAtlas, Test
using QAtlas: Bound, BB84KeyRate, Infinite

@testset "Bound(:QuantumInformation) BB84KeyRate — Shor-Preskill R = 1 - 2 H₂(e)" begin
    qi = Bound(:QuantumInformation)
    h2(e) = -e * log2(e) - (1 - e) * log2(1 - e)

    @test QAtlas.fetch(qi, BB84KeyRate(), Infinite(); qber=0.0) ≈ 1.0 atol = 1e-12
    @test QAtlas.fetch(qi, BB84KeyRate(), Infinite(); qber=0.05) ≈ 1 - 2 * h2(0.05) atol =
        1e-12
    @test_throws ArgumentError QAtlas.fetch(qi, BB84KeyRate(), Infinite(); qber=-0.1)
    @test_throws ArgumentError QAtlas.fetch(qi, BB84KeyRate(), Infinite(); qber=0.7)

    @testset "positive below, ~zero at, negative above the ~11% threshold" begin
        @test QAtlas.fetch(qi, BB84KeyRate(), Infinite(); qber=0.05) > 0
        @test QAtlas.fetch(qi, BB84KeyRate(), Infinite(); qber=0.11) ≈ 0 atol = 0.02
        @test QAtlas.fetch(qi, BB84KeyRate(), Infinite(); qber=0.2) < 0
    end

    @testset "lower bound (:geq) on the achievable key fraction" begin
        e = 0.05
        rate = 1 - 2 * h2(e)
        s = verify_bound(
            qi,
            BB84KeyRate(),
            Infinite();
            route=:saturating_constant,
            measured=[rate],
            relation=:geq,
            saturating=true,
            refs=["ShorPreskill2000"],
            fetch_kw=(; qber=e),
        )
        @test s ≈ rate atol = 1e-12
    end

    @testset "registered as a lower :bound" begin
        d = only(QAtlas.definitions(qi, BB84KeyRate(), Infinite()))
        @test d.status === :bound
        @test d.direction === :lower
    end
end
