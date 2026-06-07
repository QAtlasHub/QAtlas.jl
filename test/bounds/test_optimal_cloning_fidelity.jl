# Universal quantum cloning fidelity bound — Bound{:QuantumInformation},
# extracted from the former Universality(:QuantumMechanics) dump.
# (verify_bound is in scope via the suite's global include of test/util/verify.jl.)

using QAtlas, Test
using QAtlas: Bound, OptimalCloningFidelity, Infinite

@testset "Bound(:QuantumInformation) OptimalCloningFidelity — Buzek-Hillery F ≤ 5/6" begin
    qi = Bound(:QuantumInformation)

    @test QAtlas.fetch(qi, OptimalCloningFidelity(), Infinite()) ≈ 5 / 6 atol = 1e-14

    @testset "saturated by the optimal universal cloner" begin
        s = verify_bound(
            qi,
            OptimalCloningFidelity(),
            Infinite();
            route=:saturating_constant,
            measured=[5 / 6],
            relation=:leq,
            saturating=true,
            refs=["BuzekHillery1996"],
        )
        @test s ≈ 5 / 6 atol = 1e-14
    end

    @testset "registered as an upper :bound" begin
        d = only(QAtlas.definitions(qi, OptimalCloningFidelity(), Infinite()))
        @test d.status === :bound
        @test d.direction === :upper
        @test d.canonical
    end
end
