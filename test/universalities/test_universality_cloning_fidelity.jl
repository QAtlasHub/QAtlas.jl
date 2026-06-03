using Test
using QAtlas

@testset "Universality OptimalCloningFidelity (Buzek-Hillery 1996)" begin
    u = QAtlas.Universality(:Ising)

    @testset "Canonical qubit (d = 2): F = 5/6" begin
        f = QAtlas.fetch(u, QAtlas.OptimalCloningFidelity(), QAtlas.Infinite())
        @test f ≈ 5 / 6 atol = 1e-12
        f_explicit = QAtlas.fetch(
            u, QAtlas.OptimalCloningFidelity(), QAtlas.Infinite(); d=2
        )
        @test f_explicit ≈ 5 / 6 atol = 1e-12
    end

    @testset "Higher-dim explicit values" begin
        @test QAtlas.fetch(u, QAtlas.OptimalCloningFidelity(), QAtlas.Infinite(); d=3) ≈
            3 / 4 atol = 1e-12
        @test QAtlas.fetch(u, QAtlas.OptimalCloningFidelity(), QAtlas.Infinite(); d=4) ≈
            7 / 10 atol = 1e-12
        @test QAtlas.fetch(u, QAtlas.OptimalCloningFidelity(), QAtlas.Infinite(); d=10) ≈
            1 / 2 + 1 / 11 atol = 1e-12
    end

    @testset "Asymptotic d → ∞ gives 1/2 (classical limit)" begin
        f_large = QAtlas.fetch(
            u, QAtlas.OptimalCloningFidelity(), QAtlas.Infinite(); d=10_000
        )
        @test abs(f_large - 0.5) < 1e-3
        f_huge = QAtlas.fetch(u, QAtlas.OptimalCloningFidelity(), QAtlas.Infinite(); d=10^6)
        @test abs(f_huge - 0.5) < 1e-5
    end

    @testset "Monotonically decreasing in d" begin
        prev = 1.0
        for d in 2:20
            f = QAtlas.fetch(u, QAtlas.OptimalCloningFidelity(), QAtlas.Infinite(); d=d)
            @test f < prev
            prev = f
        end
    end

    @testset "Class-independent" begin
        vals = [
            QAtlas.fetch(
                QAtlas.Universality(sym), QAtlas.OptimalCloningFidelity(), QAtlas.Infinite()
            ) for sym in (:Ising, :Heisenberg, :XY, :Potts3, :Potts4)
        ]
        @test all(isapprox(v, 5 / 6; atol=1e-12) for v in vals)
    end

    @testset "DomainError on d < 2" begin
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.OptimalCloningFidelity(), QAtlas.Infinite(); d=1
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.OptimalCloningFidelity(), QAtlas.Infinite(); d=0
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.OptimalCloningFidelity(), QAtlas.Infinite(); d=-5
        )
    end
end
