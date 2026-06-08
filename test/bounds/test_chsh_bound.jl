# CHSH / Bell-inequality bounds — a *bound* in the bounds/ namespace
# (Bound{:QuantumInformation}), NOT a universality class.  Three theory-regime
# definitions selected by scheme=, each with its own references (whose bound).
# (verify_bound is in scope via the suite's global include of test/util/verify.jl.)

using QAtlas, Test
using QAtlas: Bound, CHSHBound, Infinite

@testset "Bound(:QuantumInformation) CHSHBound (classical 2 / quantum 2√2 / no-signalling 4)" begin
    qi = Bound(:QuantumInformation)

    @testset "scheme selects the theory regime" begin
        @test QAtlas.fetch(qi, CHSHBound(), Infinite(); scheme=:classical) ≈ 2.0 atol =
            1e-14
        @test QAtlas.fetch(qi, CHSHBound(), Infinite(); scheme=:quantum) ≈ 2sqrt(2) atol =
            1e-14
        @test QAtlas.fetch(qi, CHSHBound(), Infinite(); scheme=:no_signalling) ≈ 4.0 atol =
            1e-14
        # bare fetch returns the canonical (quantum) bound
        @test QAtlas.fetch(qi, CHSHBound(), Infinite()) ≈ 2sqrt(2) atol = 1e-14
    end

    @testset "hierarchy classical < quantum < no-signalling" begin
        b(s) = QAtlas.fetch(qi, CHSHBound(), Infinite(); scheme=s)
        @test b(:classical) < b(:quantum) < b(:no_signalling)
        @test b(:quantum) / b(:classical) ≈ sqrt(2) atol = 1e-14   # Bell-violation factor
    end

    @testset "argument validation" begin
        @test_throws ArgumentError QAtlas.fetch(qi, CHSHBound(), Infinite(); scheme=:bogus)
        # the removed `source=` selector is rejected loudly, not silently swallowed
        # (regression: it used to fall through to the default scheme and return 2√2)
        @test_throws ArgumentError QAtlas.fetch(qi, CHSHBound(), Infinite(); source=:bell)
    end

    @testset "verify_bound — quantum bound saturated by the optimal Bell state" begin
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
    end

    @testset "definitions: three :bound schemes, per-scheme refs, canonical = quantum" begin
        defs = QAtlas.definitions(qi, CHSHBound(), Infinite())
        @test length(defs) == 3
        @test all(d -> d.status === :bound && d.direction === :upper, defs)
        @test Set(d.scheme for d in defs) == Set([:classical, :quantum, :no_signalling])
        @test only(d for d in defs if d.canonical).scheme === :quantum
        @test QAtlas.validity(qi, CHSHBound(); scheme=:quantum).references ==
            ["Tsirelson1980"]
        @test QAtlas.validity(qi, CHSHBound(); scheme=:classical).references == ["CHSH1969"]
        @test QAtlas.validity(qi, CHSHBound(); scheme=:no_signalling).references ==
            ["PopescuRohrlich1994"]
    end
end
