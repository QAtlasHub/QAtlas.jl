# Mermin 3-party Bell bound — a 2nd bound in the QuantumInformation domain,
# extracted from the former Universality(:QuantumMechanics) dumping ground.
# (verify_bound is in scope via the suite's global include of test/util/verify.jl.)

using QAtlas, Test
using QAtlas: Bound, MerminGHZBound, Infinite

@testset "Bound(:QuantumInformation) MerminGHZBound (classical 2 / Mermin 4)" begin
    qi = Bound(:QuantumInformation)

    @test QAtlas.fetch(qi, MerminGHZBound(), Infinite(); source=:classical) ≈ 2.0 atol =
        1e-14
    @test QAtlas.fetch(qi, MerminGHZBound(), Infinite(); source=:mermin) ≈ 4.0 atol = 1e-14
    @test QAtlas.fetch(qi, MerminGHZBound(), Infinite()) ≈ 4.0 atol = 1e-14   # default = quantum
    @test_throws ArgumentError QAtlas.fetch(qi, MerminGHZBound(), Infinite(); source=:bogus)

    @testset "GHZ state saturates the quantum bound" begin
        ghz_value = 4.0    # |<M3>| for the GHZ state — independent witness
        s = verify_bound(
            qi,
            MerminGHZBound(),
            Infinite();
            route=:saturating_constant,
            measured=[ghz_value],
            relation=:leq,
            saturating=true,
            refs=["Mermin1990"],
        )
        @test s ≈ 4.0 atol = 1e-14
    end

    @testset "registered as an upper :bound" begin
        impl = only(
            e for e in QAtlas.REGISTRY if
            e.model === Bound{:QuantumInformation} && e.quantity === MerminGHZBound
        )
        @test impl.status === :bound
        @test impl.direction === :upper
    end
end
