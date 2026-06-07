# Sekino-Susskind fast-scrambling time — a *lower* bound (Bound{:Dynamics}),
# resolved from the deferred Universality(:?) dump.  t_* = (β/2π) log N.
# (verify_bound is in scope via the suite's global include of test/util/verify.jl.)

using QAtlas, Test
using QAtlas: Bound, ScramblingTime, Infinite

@testset "Bound(:Dynamics) ScramblingTime — Sekino-Susskind t_* = (β/2π) log N" begin
    dyn = Bound(:Dynamics)

    @test QAtlas.fetch(dyn, ScramblingTime(), Infinite(); β=2π, N=ℯ) ≈ 1.0 atol = 1e-12
    @test QAtlas.fetch(dyn, ScramblingTime(), Infinite(); β=1.0, N=100.0) ≈
        log(100.0) / (2π) atol = 1e-12
    @test_throws ArgumentError QAtlas.fetch(
        dyn, ScramblingTime(), Infinite(); β=-1.0, N=10.0
    )
    @test_throws ArgumentError QAtlas.fetch(dyn, ScramblingTime(), Infinite(); β=1.0, N=1.0)

    @testset "monotone increasing in N and β" begin
        f(β, N) = QAtlas.fetch(dyn, ScramblingTime(), Infinite(); β=β, N=N)
        @test f(1.0, 1000.0) > f(1.0, 10.0)
        @test f(2.0, 10.0) > f(1.0, 10.0)
    end

    @testset "lower bound (:geq), saturated by black holes" begin
        β, N = 1.0, 100.0
        tstar = β * log(N) / (2π)
        s = verify_bound(
            dyn,
            ScramblingTime(),
            Infinite();
            route=:saturating_constant,
            measured=[tstar],
            relation=:geq,
            saturating=true,
            refs=["SekinoSusskind2008"],
            fetch_kw=(; β=β, N=N),
        )
        @test s ≈ tstar atol = 1e-12
    end

    @testset "registered as a lower :bound" begin
        d = only(QAtlas.definitions(dyn, ScramblingTime(), Infinite()))
        @test d.status === :bound
        @test d.direction === :lower
    end
end
