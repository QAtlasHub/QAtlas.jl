# Mermin 3-party Bell bound — a 2nd quantity in the QuantumInformation domain,
# extracted from the former Universality(:QuantumMechanics) dump.  Two
# theory-regime definitions selected by scheme=.
# (verify_bound is in scope via the suite's global include of test/util/verify.jl.)

using QAtlas, Test
using QAtlas: Bound, MerminGHZBound, Infinite

@testset "Bound(:QuantumInformation) MerminGHZBound (classical 2 / quantum 4)" begin
    qi = Bound(:QuantumInformation)

    @test QAtlas.fetch(qi, MerminGHZBound(), Infinite(); scheme=:classical) ≈ 2.0 atol =
        1e-14
    @test QAtlas.fetch(qi, MerminGHZBound(), Infinite(); scheme=:quantum) ≈ 4.0 atol = 1e-14
    @test QAtlas.fetch(qi, MerminGHZBound(), Infinite()) ≈ 4.0 atol = 1e-14   # default = quantum
    @test_throws ArgumentError QAtlas.fetch(qi, MerminGHZBound(), Infinite(); scheme=:bogus)
    # the removed `source=` selector is rejected loudly, not silently swallowed
    @test_throws ArgumentError QAtlas.fetch(qi, MerminGHZBound(), Infinite(); source=:bell)

    @testset "GHZ state saturates the quantum bound" begin
        s = verify_bound(
            qi,
            MerminGHZBound(),
            Infinite();
            route=:saturating_constant,
            measured=[4.0],
            relation=:leq,
            saturating=true,
            refs=["Mermin1990"],
        )
        @test s ≈ 4.0 atol = 1e-14
    end

    @testset "definitions: two :bound schemes, canonical = quantum" begin
        defs = QAtlas.definitions(qi, MerminGHZBound(), Infinite())
        @test length(defs) == 2
        @test all(d -> d.status === :bound && d.direction === :upper, defs)
        @test only(d for d in defs if d.canonical).scheme === :quantum
    end
end
