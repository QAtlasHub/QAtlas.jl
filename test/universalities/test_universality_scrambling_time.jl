using Test
using QAtlas

@testset "Universality ScramblingTime (Sekino-Susskind 2008)" begin
    @testset "Universal formula t = (β/2π) log N" begin
        for sym in (:Ising, :Heisenberg, :XY, :Potts3, :Potts4)
            u = QAtlas.Universality(sym)
            for β in (0.5, 1.0, 5.0)
                for N in (2.0, 10.0, 100.0, 1_000_000.0)
                    t = QAtlas.fetch(
                        u, QAtlas.ScramblingTime(), QAtlas.Infinite(); beta=β, N=N
                    )
                    @test t ≈ β * log(N) / (2π) atol = 1e-12
                end
            end
        end
    end

    @testset "Class-independent at fixed (β, N)" begin
        β, N = 3.0, 100.0
        vals = [
            QAtlas.fetch(
                QAtlas.Universality(sym),
                QAtlas.ScramblingTime(),
                QAtlas.Infinite();
                beta=β,
                N=N,
            ) for sym in (:Ising, :Heisenberg, :XY)
        ]
        @test all(isapprox(v, first(vals); atol=1e-12) for v in vals)
    end

    @testset "Scaling: doubling β doubles t; squaring N doubles t" begin
        u = QAtlas.Universality(:Heisenberg)
        t_b1_N10 = QAtlas.fetch(
            u, QAtlas.ScramblingTime(), QAtlas.Infinite(); beta=1.0, N=10.0
        )
        t_b2_N10 = QAtlas.fetch(
            u, QAtlas.ScramblingTime(), QAtlas.Infinite(); beta=2.0, N=10.0
        )
        t_b1_N100 = QAtlas.fetch(
            u, QAtlas.ScramblingTime(), QAtlas.Infinite(); beta=1.0, N=100.0
        )
        @test t_b2_N10 ≈ 2 * t_b1_N10 atol = 1e-12
        @test t_b1_N100 ≈ 2 * t_b1_N10 atol = 1e-12
    end

    @testset "DomainError on invalid args" begin
        u = QAtlas.Universality(:Ising)
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.ScramblingTime(), QAtlas.Infinite(); beta=0.0, N=10.0
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.ScramblingTime(), QAtlas.Infinite(); beta=-1.0, N=10.0
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.ScramblingTime(), QAtlas.Infinite(); beta=1.0, N=1.0
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.ScramblingTime(), QAtlas.Infinite(); beta=1.0, N=0.5
        )
    end
end
