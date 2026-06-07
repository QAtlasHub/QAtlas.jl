using Test
using QAtlas
using QAtlas:
    Universality, LogarithmicNegativity, MutualInformation, Infinite, CentralCharge

@testset "Universality LogarithmicNegativity adjacent intervals (CCT 2012, #580)" begin
    @testset "Closed form (c/4) log[ℓ_A ℓ_B / (ℓ_A + ℓ_B)] across 5 classes" begin
        d_for = Dict(:Ising=>2, :XY=>2, :Heisenberg=>1, :Potts3=>2, :Potts4=>2)
        for class in (:Ising, :XY, :Heisenberg, :Potts3, :Potts4)
            c = float(QAtlas.fetch(Universality(class), CentralCharge(); d=d_for[class]))
            for (ℓ_A, ℓ_B) in ((5.0, 10.0), (1.0, 100.0), (50.0, 50.0))
                E = QAtlas.fetch(
                    Universality(class),
                    LogarithmicNegativity(),
                    Infinite();
                    ℓ_A=ℓ_A,
                    ℓ_B=ℓ_B,
                )
                expected = (c / 4) * log(ℓ_A * ℓ_B / (ℓ_A + ℓ_B))
                @test E ≈ expected atol = 1e-12
            end
        end
    end

    @testset "Ratio E_neg / I_MI = 3/4 universally" begin
        # The negativity differs from mutual info only by prefactor c/4
        # vs c/3, so their ratio is exactly 3/4 across all classes and
        # all (ℓ_A, ℓ_B) values.
        for class in (:Ising, :XY, :Heisenberg, :Potts3, :Potts4)
            for (ℓ_A, ℓ_B) in ((5.0, 10.0), (20.0, 50.0), (100.0, 200.0))
                E = QAtlas.fetch(
                    Universality(class),
                    LogarithmicNegativity(),
                    Infinite();
                    ℓ_A=ℓ_A,
                    ℓ_B=ℓ_B,
                )
                I = QAtlas.fetch(
                    Universality(class), MutualInformation(), Infinite(); ℓ_A=ℓ_A, ℓ_B=ℓ_B
                )
                @test E / I ≈ 3 / 4 atol = 1e-12
            end
        end
    end

    @testset "E_neg > 0 for ℓ_A · ℓ_B / (ℓ_A + ℓ_B) > 1" begin
        for class in (:Ising, :Heisenberg)
            for (ℓ_A, ℓ_B) in ((5.0, 10.0), (50.0, 100.0))
                E = QAtlas.fetch(
                    Universality(class),
                    LogarithmicNegativity(),
                    Infinite();
                    ℓ_A=ℓ_A,
                    ℓ_B=ℓ_B,
                )
                @test E > 0
            end
        end
    end

    @testset "Argument validation" begin
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), LogarithmicNegativity(), Infinite(); ℓ_A=0.0, ℓ_B=5.0
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), LogarithmicNegativity(), Infinite(); ℓ_A=5.0, ℓ_B=-1.0
        )
    end
end
