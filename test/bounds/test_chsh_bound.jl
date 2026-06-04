# CHSH / Bell-inequality bounds — proof that a *bound* lives in the bounds/
# namespace (Bound{:QuantumInformation}), NOT under a universality class.
# Exercises: the source= selector (whose bound), the classical < quantum <
# no-signalling ladder, verify_bound saturation, and the registry
# direction=:upper metadata.  (verify_bound is in scope via the suite's
# global include of test/util/verify.jl, like test_TFIM_status_examples.jl.)

using QAtlas, Test
using QAtlas: Bound, CHSHBound, Infinite

@testset "Bound(:QuantumInformation) CHSHBound (Bell 2 / Tsirelson 2√2 / PR 4)" begin
    qi = Bound(:QuantumInformation)

    @testset "source selects whose bound" begin
        @test QAtlas.fetch(qi, CHSHBound(), Infinite(); source=:bell) ≈ 2.0 atol = 1e-14
        @test QAtlas.fetch(qi, CHSHBound(), Infinite(); source=:tsirelson) ≈ 2sqrt(2) atol =
            1e-14
        @test QAtlas.fetch(qi, CHSHBound(), Infinite(); source=:popescu_rohrlich) ≈ 4.0 atol =
            1e-14
        # default is the quantum (Tsirelson) bound
        @test QAtlas.fetch(qi, CHSHBound(), Infinite()) ≈ 2sqrt(2) atol = 1e-14
    end

    @testset "hierarchy classical < quantum < no-signalling" begin
        b(src) = QAtlas.fetch(qi, CHSHBound(), Infinite(); source=src)
        @test b(:bell) < b(:tsirelson) < b(:popescu_rohrlich)
        @test b(:tsirelson) / b(:bell) ≈ sqrt(2) atol = 1e-14   # Bell-violation factor
    end

    @testset "argument validation" begin
        @test_throws ArgumentError QAtlas.fetch(qi, CHSHBound(), Infinite(); source=:bogus)
    end

    @testset "verify_bound — the optimal Bell state saturates the quantum bound" begin
        optimal = 2sqrt(2)            # CHSH of the optimal Bell state — independent witness
        s = verify_bound(
            qi,
            CHSHBound(),
            Infinite();
            route=:saturating_constant,
            measured=[optimal],
            relation=:leq,
            saturating=true,
            refs=["Tsirelson1980"],
        )
        @test s ≈ 2sqrt(2) atol = 1e-14
        @test optimal <= s + 1e-12
    end

    @testset "registered as an upper :bound, not a universality class" begin
        impls = [
            e for e in QAtlas.REGISTRY if
            e.model === Bound{:QuantumInformation} && e.quantity === CHSHBound
        ]
        @test length(impls) == 1
        @test impls[1].status === :bound
        @test impls[1].direction === :upper       # pinned in the registry (FW-1)
    end
end
