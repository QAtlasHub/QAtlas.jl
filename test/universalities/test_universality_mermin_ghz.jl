using Test
using QAtlas

@testset "Universality MerminGHZBound (Mermin 1990)" begin
    @testset "Classical bound = 2" begin
        for sym in (:Ising, :Heisenberg, :XY, :Potts3, :Potts4)
            u = QAtlas.Universality(sym)
            @test QAtlas.fetch(
                u, QAtlas.MerminGHZBound(), QAtlas.Infinite(); theory=:classical
            ) == 2.0
        end
    end

    @testset "Quantum bound = 4 (GHZ saturates)" begin
        for sym in (:Ising, :Heisenberg, :XY)
            u = QAtlas.Universality(sym)
            @test QAtlas.fetch(
                u, QAtlas.MerminGHZBound(), QAtlas.Infinite(); theory=:quantum
            ) == 4.0
        end
    end

    @testset "No-signalling bound = 4 (QM saturates)" begin
        u = QAtlas.Universality(:Ising)
        @test QAtlas.fetch(
            u, QAtlas.MerminGHZBound(), QAtlas.Infinite(); theory=:no_signalling
        ) == 4.0
    end

    @testset "QM bound = no-signalling bound (unlike CHSH)" begin
        u = QAtlas.Universality(:Heisenberg)
        @test QAtlas.fetch(
            u, QAtlas.MerminGHZBound(), QAtlas.Infinite(); theory=:quantum
        ) == QAtlas.fetch(
            u, QAtlas.MerminGHZBound(), QAtlas.Infinite(); theory=:no_signalling
        )
    end

    @testset "Classical / quantum ratio = 1/2 (factor of 2 gap)" begin
        u = QAtlas.Universality(:Ising)
        c = QAtlas.fetch(u, QAtlas.MerminGHZBound(), QAtlas.Infinite(); theory=:classical)
        q = QAtlas.fetch(u, QAtlas.MerminGHZBound(), QAtlas.Infinite(); theory=:quantum)
        @test c / q ≈ 0.5 atol = 1e-12
    end

    @testset "Universality-class independent" begin
        vals_q = [
            QAtlas.fetch(
                QAtlas.Universality(sym),
                QAtlas.MerminGHZBound(),
                QAtlas.Infinite();
                theory=:quantum,
            ) for sym in (:Ising, :Heisenberg, :XY, :Potts3, :Potts4)
        ]
        @test all(v == 4.0 for v in vals_q)
    end

    @testset "ArgumentError for unknown theory" begin
        u = QAtlas.Universality(:Ising)
        @test_throws ArgumentError QAtlas.fetch(
            u, QAtlas.MerminGHZBound(), QAtlas.Infinite(); theory=:nonsense
        )
        @test_throws ArgumentError QAtlas.fetch(
            u, QAtlas.MerminGHZBound(), QAtlas.Infinite(); theory=:hidden
        )
    end
end
