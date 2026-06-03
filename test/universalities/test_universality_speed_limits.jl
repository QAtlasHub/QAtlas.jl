using Test
using QAtlas

@testset "Universality QuantumSpeedLimits (Mandelstam-Tamm + Margolus-Levitin)" begin
    @testset "MandelstamTammBound = π / (2 ΔE)" begin
        for sym in (:Ising, :Heisenberg, :XY, :Potts3, :Potts4)
            u = QAtlas.Universality(sym)
            for ΔE in (0.5, 1.0, 2.5, 10.0, 100.0)
                t = QAtlas.fetch(
                    u, QAtlas.MandelstamTammBound(), QAtlas.Infinite(); delta_E=ΔE
                )
                @test t ≈ π / (2 * ΔE) atol = 1e-12
            end
        end
    end

    @testset "MargolusLevitinBound = π / (2 E)" begin
        for sym in (:Ising, :Heisenberg, :XY, :Potts3, :Potts4)
            u = QAtlas.Universality(sym)
            for E in (0.5, 1.0, 2.5, 10.0, 100.0)
                t = QAtlas.fetch(
                    u, QAtlas.MargolusLevitinBound(), QAtlas.Infinite(); mean_E=E
                )
                @test t ≈ π / (2 * E) atol = 1e-12
            end
        end
    end

    @testset "MT = ML when ΔE = E - E_0 (coincidence point)" begin
        u = QAtlas.Universality(:Heisenberg)
        for x in (0.7, 1.3, 5.0)
            t_MT = QAtlas.fetch(
                u, QAtlas.MandelstamTammBound(), QAtlas.Infinite(); delta_E=x
            )
            t_ML = QAtlas.fetch(
                u, QAtlas.MargolusLevitinBound(), QAtlas.Infinite(); mean_E=x
            )
            @test t_MT ≈ t_ML atol = 1e-12
        end
    end

    @testset "Class independence" begin
        ΔE = 3.0
        E = 2.0
        for sym in (:Ising, :Heisenberg, :XY)
            u = QAtlas.Universality(sym)
            t_MT = QAtlas.fetch(
                u, QAtlas.MandelstamTammBound(), QAtlas.Infinite(); delta_E=ΔE
            )
            t_ML = QAtlas.fetch(
                u, QAtlas.MargolusLevitinBound(), QAtlas.Infinite(); mean_E=E
            )
            @test t_MT ≈ π / 6 atol = 1e-12
            @test t_ML ≈ π / 4 atol = 1e-12
        end
    end

    @testset "DomainError on non-positive energy" begin
        u = QAtlas.Universality(:Ising)
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.MandelstamTammBound(), QAtlas.Infinite(); delta_E=0.0
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.MandelstamTammBound(), QAtlas.Infinite(); delta_E=-1.0
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.MargolusLevitinBound(), QAtlas.Infinite(); mean_E=0.0
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.MargolusLevitinBound(), QAtlas.Infinite(); mean_E=-2.0
        )
    end
end
