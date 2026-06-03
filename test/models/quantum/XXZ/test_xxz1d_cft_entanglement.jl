using Test
using QAtlas
using QAtlas: XXZ1D, VonNeumannEntropy, RenyiEntropy, Universality, Infinite

@testset "XXZ1D Infinite VN/Renyi delegation to Universality (#580 Phase 2)" begin
    @testset "Critical regime -1 < Δ < 1 delegates to Universality(:XY) (c=1)" begin
        for Δ in (-0.9, -0.5, -0.1, 0.0, 0.3, 0.7, 0.99)
            m = XXZ1D(; J=1.0, Δ=Δ)
            for ℓ in (5, 20, 100), β in (Inf, 50.0, 5.0)
                S_model = QAtlas.fetch(m, VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=β)
                S_univ = QAtlas.fetch(
                    Universality(:XY), VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=β
                )
                @test S_model ≈ S_univ atol = 1e-14
            end
        end
    end

    @testset "Δ = 1 SU(2) Heisenberg routes via Universality(:Heisenberg) (c=1)" begin
        m = XXZ1D(; J=1.0, Δ=1.0)
        for ℓ in (10, 50), β in (Inf, 10.0, 1.0)
            S_model = QAtlas.fetch(m, VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=β)
            S_univ = QAtlas.fetch(
                Universality(:Heisenberg), VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=β
            )
            @test S_model ≈ S_univ atol = 1e-14
        end
    end

    @testset "Renyi-α critical regime" begin
        for Δ in (-0.5, 0.0, 0.5, 1.0), α in (0.5, 2.0, 5.0)
            m = XXZ1D(; J=1.0, Δ=Δ)
            for ℓ in (10, 30), β in (Inf, 5.0)
                S = QAtlas.fetch(m, RenyiEntropy(α), Infinite(); ℓ=ℓ, beta=β)
                class = Δ == 1.0 ? :Heisenberg : :XY
                S_u = QAtlas.fetch(
                    Universality(class), RenyiEntropy(α), Infinite(); ℓ=ℓ, beta=β
                )
                @test S ≈ S_u atol = 1e-14
            end
        end
    end

    @testset "T = 0 critical c = 1 leading log: (1/3) log ℓ across critical regime" begin
        for Δ in (-0.9, -0.5, 0.0, 0.5, 1.0)
            m = XXZ1D(; J=1.0, Δ=Δ)
            for ℓ in (5, 20, 50)
                S = QAtlas.fetch(m, VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=Inf)
                @test S ≈ (1 / 3) * log(ℓ) atol = 1e-12
            end
        end
    end

    @testset "Gapped |Δ| > 1 throws DomainError" begin
        for Δ in (-2.0, -1.5, 1.5, 2.0, 5.0)
            m = XXZ1D(; J=1.0, Δ=Δ)
            @test_throws DomainError QAtlas.fetch(
                m, VonNeumannEntropy(), Infinite(); ℓ=10, beta=Inf
            )
            @test_throws DomainError QAtlas.fetch(
                m, RenyiEntropy(2.0), Infinite(); ℓ=10, beta=Inf
            )
        end
    end
end
