using Test
using QAtlas
using QAtlas: Universality, PageEntropy, Infinite

@testset "Universality(:HaarRandom) PageEntropy Page 1993 (#580)" begin
    @testset "m = 1: subsystem A is a single state -> S = 0" begin
        for n in (1, 2, 4, 10, 100)
            S = QAtlas.fetch(
                Universality(:HaarRandom), PageEntropy(), Infinite(); d_A=1, d_B=n
            )
            @test S ≈ 0.0 atol = 1e-14
        end
    end

    @testset "m = n = 2: <S_A> = sum_{k=3}^{4} 1/k - 1/4 = 7/12" begin
        S = QAtlas.fetch(Universality(:HaarRandom), PageEntropy(), Infinite(); d_A=2, d_B=2)
        # 1/3 + 1/4 - 1/4 = 1/3
        @test S ≈ 1 / 3 atol = 1e-14
    end

    @testset "A ↔ B symmetry (purity of global state)" begin
        for (dA, dB) in ((2, 4), (3, 8), (4, 16), (5, 7))
            S1 = QAtlas.fetch(
                Universality(:HaarRandom), PageEntropy(), Infinite(); d_A=dA, d_B=dB
            )
            S2 = QAtlas.fetch(
                Universality(:HaarRandom), PageEntropy(), Infinite(); d_A=dB, d_B=dA
            )
            @test S1 ≈ S2 atol = 1e-14
        end
    end

    @testset "m = n asymptotic: <S_A> → log m - 1/2 (large m)" begin
        # Confirm convergence to Page asymptotic at successively larger m.
        for m in (8, 16, 32, 64)
            S = QAtlas.fetch(
                Universality(:HaarRandom), PageEntropy(), Infinite(); d_A=m, d_B=m
            )
            @test abs(S - (log(m) - 0.5)) < 1.0 / m
        end
    end

    @testset "Sub-maximal: 0 ≤ <S_A> ≤ log(min(d_A, d_B))" begin
        for dA in (1, 2, 4, 8), dB in (1, 2, 4, 8, 16)
            S = QAtlas.fetch(
                Universality(:HaarRandom), PageEntropy(), Infinite(); d_A=dA, d_B=dB
            )
            @test S >= 0
            @test S <= log(min(dA, dB)) + 1e-14
        end
    end

    @testset "Argument validation" begin
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:HaarRandom), PageEntropy(), Infinite(); d_A=0, d_B=4
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:HaarRandom), PageEntropy(), Infinite(); d_A=4, d_B=-1
        )
    end
end
