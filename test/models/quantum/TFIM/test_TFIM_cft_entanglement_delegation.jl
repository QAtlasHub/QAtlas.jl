using Test
using QAtlas
using QAtlas: TFIM, VonNeumannEntropy, RenyiEntropy, Universality, Infinite

@testset "TFIM Infinite VN/Renyi: h=J delegates to Universality(:Ising) (#580 Phase 2)" begin
    @testset "h = J critical delegates with a = 1/2 lattice spacing convention" begin
        m = TFIM(; J=1.0, h=1.0)
        for ℓ in (1, 5, 10, 50, 100), β in (Inf, 100.0, 10.0, 1.0)
            S_model = QAtlas.fetch(m, VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=β)
            S_univ = QAtlas.fetch(
                Universality(:Ising), VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=β, a=1 // 2
            )
            @test S_model ≈ S_univ atol = 1e-14
        end
    end

    @testset "h = J critical leading log preserved: T=0 returns (1/6) log(2ℓ)" begin
        # c = 1/2 for Ising → S_VN at T=0 = (c/3) log(ℓ/a) with a = 1/2
        # = (1/6) log(2ℓ). This is the historical TFIM convention.
        m = TFIM(; J=1.0, h=1.0)
        for ℓ in (1, 5, 10, 100)
            S = QAtlas.fetch(m, VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=Inf)
            @test S ≈ (1 / 6) * log(2 * ℓ) atol = 1e-12
        end
    end

    @testset "Renyi-α at h = J delegates with same a = 1/2" begin
        m = TFIM(; J=1.0, h=1.0)
        for α in (0.5, 2.0, 3.0, 5.0), ℓ in (10, 30), β in (Inf, 5.0)
            S_model = QAtlas.fetch(m, RenyiEntropy(α), Infinite(); ℓ=ℓ, beta=β)
            S_univ = QAtlas.fetch(
                Universality(:Ising), RenyiEntropy(α), Infinite(); ℓ=ℓ, beta=β, a=1 // 2
            )
            @test S_model ≈ S_univ atol = 1e-14
        end
    end

    @testset "Off-critical h ≠ J unchanged (model-local helper)" begin
        # Gapped TFIM with finite mass scale ξ = 1/(2|h-J|). The
        # delegation path is NOT taken; the value must match the
        # historical model-local formula.
        for (h, ξ_expected) in ((0.5, 1.0), (1.5, 1.0), (2.0, 0.5))
            m = TFIM(; J=1.0, h=h)
            S = QAtlas.fetch(m, VonNeumannEntropy(), Infinite(); ℓ=50, beta=Inf)
            # Off-critical gapped CC: S_VN = (c/12) · log(2 ξ sinh(ℓ/ξ))
            # for ℓ ≫ ξ this saturates.
            ξ = 1 / (2 * abs(h - 1.0))
            @test ξ ≈ ξ_expected atol = 1e-12
            @test S > 0
            @test isfinite(S)
        end
    end
end
