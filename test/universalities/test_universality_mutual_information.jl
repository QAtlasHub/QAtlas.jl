using Test
using QAtlas
using QAtlas: Universality, MutualInformation, VonNeumannEntropy, Infinite, CentralCharge

@testset "Universality MutualInformation adjacent intervals (#580)" begin
    @testset "T=0 closed form: I = (c/3) log[ℓ_A ℓ_B / (ℓ_A + ℓ_B)]" begin
        d_for = Dict(:Ising=>2, :XY=>2, :Heisenberg=>1, :Potts3=>2, :Potts4=>2)
        for class in (:Ising, :XY, :Heisenberg, :Potts3, :Potts4)
            c = float(QAtlas.fetch(Universality(class), CentralCharge(); d=d_for[class]))
            for (ℓ_A, ℓ_B) in ((5.0, 10.0), (1.0, 1.0), (50.0, 100.0))
                I = QAtlas.fetch(
                    Universality(class), MutualInformation(), Infinite(); ℓ_A=ℓ_A, ℓ_B=ℓ_B
                )
                expected = (c / 3) * log(ℓ_A * ℓ_B / (ℓ_A + ℓ_B))
                @test I ≈ expected atol = 1e-12
            end
        end
    end

    @testset "I = S(A) + S(B) - S(A∪B) consistency at T = 0 + finite β" begin
        for class in (:Ising, :Heisenberg)
            for β in (Inf, 50.0, 5.0), (ℓ_A, ℓ_B) in ((3.0, 7.0), (10.0, 20.0))
                I = QAtlas.fetch(
                    Universality(class),
                    MutualInformation(),
                    Infinite();
                    ℓ_A=ℓ_A,
                    ℓ_B=ℓ_B,
                    beta=β,
                )
                S_A = QAtlas.fetch(
                    Universality(class), VonNeumannEntropy(), Infinite(); ℓ=ℓ_A, beta=β
                )
                S_B = QAtlas.fetch(
                    Universality(class), VonNeumannEntropy(), Infinite(); ℓ=ℓ_B, beta=β
                )
                S_AB = QAtlas.fetch(
                    Universality(class),
                    VonNeumannEntropy(),
                    Infinite();
                    ℓ=ℓ_A + ℓ_B,
                    beta=β,
                )
                @test I ≈ (S_A + S_B - S_AB) atol = 1e-12
            end
        end
    end

    @testset "I > 0 for ℓ_A · ℓ_B / (ℓ_A + ℓ_B) > 1 (positive universal log)" begin
        # Only the universal log piece is returned; the non-universal
        # S_0 offsets are dropped. So I_universal > 0 iff the argument
        # of the log exceeds the UV cutoff (here a = 1). At very small
        # ℓ (ℓ_A·ℓ_B / (ℓ_A+ℓ_B) < 1) the universal piece alone can be
        # negative, with the physical I > 0 requiring the dropped
        # constant.
        for class in (:Ising, :Heisenberg),
            (ℓ_A, ℓ_B) in ((5.0, 10.0), (20.0, 50.0), (100.0, 200.0))

            I = QAtlas.fetch(
                Universality(class), MutualInformation(), Infinite(); ℓ_A=ℓ_A, ℓ_B=ℓ_B
            )
            @test I > 0
        end
    end

    @testset "Argument validation" begin
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), MutualInformation(), Infinite(); ℓ_A=-1.0, ℓ_B=5.0
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), MutualInformation(), Infinite(); ℓ_A=5.0, ℓ_B=0.0
        )
    end
end
