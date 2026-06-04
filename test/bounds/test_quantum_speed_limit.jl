# Margolus-Levitin quantum speed limit — a *lower* bound (Bound{:Dynamics}),
# extracted from the former Universality(:QuantumMechanics) dump.  Exercises the
# :lower direction (verify_bound relation=:geq).
# (verify_bound is in scope via the suite's global include of test/util/verify.jl.)

using QAtlas, Test
using QAtlas: Bound, QuantumSpeedLimit, Infinite

@testset "Bound(:Dynamics) QuantumSpeedLimit — Margolus-Levitin τ ≥ π/(2E)" begin
    dyn = Bound(:Dynamics)

    @test QAtlas.fetch(dyn, QuantumSpeedLimit(), Infinite(); E=1.0) ≈ π / 2 atol = 1e-12
    @test QAtlas.fetch(
        dyn, QuantumSpeedLimit(), Infinite(); scheme=:margolus_levitin, E=2.0
    ) ≈ π / 4 atol = 1e-12
    @test_throws ArgumentError QAtlas.fetch(dyn, QuantumSpeedLimit(), Infinite(); E=-1.0)

    @testset "lower bound (:geq), saturated by an equal two-level superposition" begin
        E = 1.0
        s = verify_bound(
            dyn,
            QuantumSpeedLimit(),
            Infinite();
            route=:saturating_constant,
            measured=[π / (2E)],
            relation=:geq,
            saturating=true,
            refs=["MargolusLevitin1998"],
            fetch_kw=(; E=E),
        )
        @test s ≈ π / 2 atol = 1e-12
    end

    @testset "registered as a lower :bound" begin
        d = only(QAtlas.definitions(dyn, QuantumSpeedLimit(), Infinite()))
        @test d.status === :bound
        @test d.direction === :lower
        @test d.canonical
    end
end
