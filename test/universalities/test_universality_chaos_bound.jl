using Test
using QAtlas

@testset "Universality ChaosBound (Maldacena-Shenker-Stanford 2016)" begin
    @testset "Bound value 2 pi / beta" begin
        for sym in (:Ising, :Heisenberg, :XY, :Potts3, :Potts4)
            u = QAtlas.Universality(sym)
            for β in (0.5, 1.0, 2.0, 5.0, 10.0)
                λ = QAtlas.fetch(u, QAtlas.ChaosBound(), QAtlas.Infinite(); beta=β)
                @test λ ≈ 2π / β atol = 1e-12
            end
        end
    end

    @testset "Universality-independent (same value across classes at fixed β)" begin
        β = 3.0
        vals = [
            QAtlas.fetch(
                QAtlas.Universality(sym), QAtlas.ChaosBound(), QAtlas.Infinite(); beta=β
            ) for sym in (:Ising, :Heisenberg, :XY, :Potts3, :Potts4)
        ]
        @test all(isapprox(v, first(vals); atol=1e-12) for v in vals)
    end

    @testset "DomainError on non-positive beta" begin
        u = QAtlas.Universality(:Ising)
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.ChaosBound(), QAtlas.Infinite(); beta=0.0
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.ChaosBound(), QAtlas.Infinite(); beta=-1.0
        )
    end

    @testset "High-T limit blows up, low-T limit vanishes" begin
        u = QAtlas.Universality(:Heisenberg)
        λ_hot = QAtlas.fetch(u, QAtlas.ChaosBound(), QAtlas.Infinite(); beta=0.01)
        λ_cold = QAtlas.fetch(u, QAtlas.ChaosBound(), QAtlas.Infinite(); beta=100.0)
        @test λ_hot > 100
        @test λ_cold < 0.1
    end
end
