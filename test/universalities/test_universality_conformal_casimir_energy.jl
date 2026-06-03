using Test
using QAtlas
using QAtlas: Universality, ConformalCasimirEnergy, Infinite, CentralCharge

@testset "Universality ConformalCasimirEnergy = -π c / (6 L) Cardy 1986 (#580)" begin
    @testset "Formula across 5 universality classes" begin
        d_for = Dict(:Ising=>2, :XY=>2, :Heisenberg=>1, :Potts3=>2, :Potts4=>2)
        for class in (:Ising, :XY, :Heisenberg, :Potts3, :Potts4)
            c = float(QAtlas.fetch(Universality(class), CentralCharge(); d=d_for[class]))
            for L in (1.0, 5.0, 20.0, 100.0, 1000.0)
                E0 = QAtlas.fetch(
                    Universality(class), ConformalCasimirEnergy(), Infinite(); L=L
                )
                expected = -π * c / (6 * L)
                @test E0 ≈ expected atol = 1e-14
            end
        end
    end

    @testset "E_0 < 0 for unitary classes (c > 0)" begin
        for class in (:Ising, :XY, :Heisenberg, :Potts3, :Potts4)
            E0 = QAtlas.fetch(
                Universality(class), ConformalCasimirEnergy(), Infinite(); L=10.0
            )
            @test E0 < 0
        end
    end

    @testset "E_0(L) decays as 1/L" begin
        for class in (:Ising, :Heisenberg)
            E1 = QAtlas.fetch(
                Universality(class), ConformalCasimirEnergy(), Infinite(); L=1.0
            )
            E10 = QAtlas.fetch(
                Universality(class), ConformalCasimirEnergy(), Infinite(); L=10.0
            )
            E100 = QAtlas.fetch(
                Universality(class), ConformalCasimirEnergy(), Infinite(); L=100.0
            )
            @test E10 ≈ E1 / 10 atol = 1e-14
            @test E100 ≈ E1 / 100 atol = 1e-14
        end
    end

    @testset "E_0 linear in c at fixed L" begin
        L = 5.0
        E_Ising = QAtlas.fetch(
            Universality(:Ising), ConformalCasimirEnergy(), Infinite(); L=L
        )
        E_Heis = QAtlas.fetch(
            Universality(:Heisenberg), ConformalCasimirEnergy(), Infinite(); L=L
        )
        c_I = float(QAtlas.fetch(Universality(:Ising), CentralCharge(); d=2))
        c_H = float(QAtlas.fetch(Universality(:Heisenberg), CentralCharge(); d=1))
        @test E_Heis / E_Ising ≈ c_H / c_I atol = 1e-12
    end

    @testset "Argument validation" begin
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), ConformalCasimirEnergy(), Infinite(); L=0.0
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Heisenberg), ConformalCasimirEnergy(), Infinite(); L=-1.0
        )
    end
end
