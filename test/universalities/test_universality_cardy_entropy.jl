using Test
using QAtlas
using QAtlas: Universality, CardyEntropy, Infinite, CentralCharge

@testset "Universality CardyEntropy (Cardy 1986, #580)" begin
    @testset "Formula 2 π sqrt(c E / 6) across 5 universality classes" begin
        d_for = Dict(:Ising=>2, :XY=>2, :Heisenberg=>1, :Potts3=>2, :Potts4=>2)
        for class in (:Ising, :XY, :Heisenberg, :Potts3, :Potts4)
            c = float(QAtlas.fetch(Universality(class), CentralCharge(); d=d_for[class]))
            for E in (0.1, 1.0, 10.0, 100.0, 1000.0)
                S = QAtlas.fetch(Universality(class), CardyEntropy(), Infinite(); E=E)
                expected = 2 * π * sqrt(c * E / 6)
                @test S ≈ expected atol = 1e-14
            end
        end
    end

    @testset "E = 0 -> S = 0 (degenerate ground state)" begin
        for class in (:Ising, :Heisenberg, :Potts3)
            S = QAtlas.fetch(Universality(class), CardyEntropy(), Infinite(); E=0.0)
            @test S ≈ 0.0 atol = 1e-14
        end
    end

    @testset "S grows as sqrt(E) (asymptotic scaling)" begin
        for class in (:Ising, :Heisenberg)
            S1 = QAtlas.fetch(Universality(class), CardyEntropy(), Infinite(); E=1.0)
            S100 = QAtlas.fetch(Universality(class), CardyEntropy(), Infinite(); E=100.0)
            @test S100 ≈ 10 * S1 atol = 1e-12
        end
    end

    @testset "S grows as sqrt(c) (universality scaling)" begin
        E = 4.0
        S_Ising = QAtlas.fetch(Universality(:Ising), CardyEntropy(), Infinite(); E=E)
        S_Heisenberg = QAtlas.fetch(
            Universality(:Heisenberg), CardyEntropy(), Infinite(); E=E
        )
        c_I = float(QAtlas.fetch(Universality(:Ising), CentralCharge(); d=2))
        c_H = float(QAtlas.fetch(Universality(:Heisenberg), CentralCharge(); d=1))
        @test S_Heisenberg / S_Ising ≈ sqrt(c_H / c_I) atol = 1e-12
    end

    @testset "Argument validation" begin
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), CardyEntropy(), Infinite(); E=-1.0
        )
    end
end
