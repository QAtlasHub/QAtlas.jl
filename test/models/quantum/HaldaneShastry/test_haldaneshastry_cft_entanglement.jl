using Test
using QAtlas
using QAtlas: HaldaneShastry, VonNeumannEntropy, RenyiEntropy, Universality, Infinite

@testset "HaldaneShastry Infinite VN/Renyi delegation to Universality(:Heisenberg) (#580 Phase 2)" begin
    m = HaldaneShastry()
    @testset "VN at (ℓ, β) matches universality output" begin
        for ℓ in (5, 10, 50, 100), β in (Inf, 100.0, 10.0, 1.0)
            S_model = QAtlas.fetch(m, VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=β)
            S_univ = QAtlas.fetch(
                Universality(:Heisenberg), VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=β
            )
            @test S_model ≈ S_univ atol = 1e-14
        end
    end

    @testset "Renyi-α at (ℓ, β) matches universality output" begin
        for α in (0.5, 2.0, 3.0, 5.0), ℓ in (10, 30), β in (Inf, 5.0)
            S_model = QAtlas.fetch(m, RenyiEntropy(α), Infinite(); ℓ=ℓ, beta=β)
            S_univ = QAtlas.fetch(
                Universality(:Heisenberg), RenyiEntropy(α), Infinite(); ℓ=ℓ, beta=β
            )
            @test S_model ≈ S_univ atol = 1e-14
        end
    end

    @testset "T=0 critical c=1 leading log: (1/3) log ℓ" begin
        for ℓ in (1, 5, 10, 50, 100)
            S = QAtlas.fetch(m, VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=Inf)
            @test S ≈ (1 / 3) * log(ℓ) atol = 1e-12
        end
    end
end
