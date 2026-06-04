# MSS chaos bound — a dynamical bound in the bounds/ namespace
# (Bound{:Dynamics}), extracted from the former Universality(:QuantumMechanics) dump.
# (verify_bound is in scope via the suite's global include of test/util/verify.jl.)

using QAtlas, Test
using QAtlas: Bound, ChaosBound, Infinite

@testset "Bound(:Dynamics) ChaosBound — MSS λ_L ≤ 2π/β" begin
    dyn = Bound(:Dynamics)

    for β in (0.5, 1.0, 4.0)
        @test QAtlas.fetch(dyn, ChaosBound(), Infinite(); β=β) ≈ 2π / β atol = 1e-12
    end
    @test_throws ArgumentError QAtlas.fetch(dyn, ChaosBound(), Infinite(); β=-1.0)

    @testset "saturated by a maximally chaotic system" begin
        β = 2.0
        s = verify_bound(
            dyn,
            ChaosBound(),
            Infinite();
            route=:saturating_constant,
            measured=[2π / β],
            relation=:leq,
            saturating=true,
            refs=["MaldacenaShenkerStanford2016"],
            fetch_kw=(; β=β),
        )
        @test s ≈ 2π / β atol = 1e-12
    end

    @testset "registered as an upper :bound" begin
        d = only(QAtlas.definitions(dyn, ChaosBound(), Infinite()))
        @test d.status === :bound
        @test d.direction === :upper
        @test d.canonical
    end
end
