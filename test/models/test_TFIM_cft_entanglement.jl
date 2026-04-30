using QAtlas, Test

@testset "TFIM CC entanglement at Infinite" begin
    @testset "Critical h = J: VN scales as (c/3) log(2 ℓ)" begin
        model = TFIM(; J=1.0, h=1.0)
        c = 0.5
        for ℓ in (10, 50, 200)
            S = QAtlas.fetch(model, VonNeumannEntropy(), Infinite(); ℓ=ℓ)
            @test S ≈ (c / 3) * log(2 * ℓ) atol = 1e-10
        end
    end

    @testset "Critical: Renyi(α) scales as (c/6)(1 + 1/α) log(2 ℓ)" begin
        model = TFIM(; J=1.0, h=1.0)
        c = 0.5
        for α in (0.5, 2.0, 3.0), ℓ in (10, 50, 200)
            S = QAtlas.fetch(model, RenyiEntropy(α), Infinite(); ℓ=ℓ)
            @test S ≈ (c / 6) * (1 + 1 / α) * log(2 * ℓ) atol = 1e-10
        end
    end

    @testset "Off-critical: gapped CC closed form" begin
        # h = 0.5, J = 1.0 → gapped, ξ = 1/(2·0.5) = 1.
        model = TFIM(; J=1.0, h=0.5)
        c = 0.5
        ξ = 1.0
        S1 = QAtlas.fetch(model, VonNeumannEntropy(), Infinite(); ℓ=1)
        S50 = QAtlas.fetch(model, VonNeumannEntropy(), Infinite(); ℓ=50)
        @test S50 > S1
        @test isfinite(S50)
        # Direct formula check (volume-like scaling once ℓ ≫ ξ; the
        # universal CC envelope still applies).
        S_predicted = (c / 6) * log(2 * ξ * sinh(50 / ξ))
        @test S50 ≈ S_predicted atol = 1e-10
    end

    @testset "ξ → ∞ continuity (h → J) at fixed ℓ" begin
        ℓ = 50
        S_cr = QAtlas.fetch(TFIM(; J=1.0, h=1.0), VonNeumannEntropy(), Infinite(); ℓ=ℓ)
        for δ in (1e-3, 1e-4, 1e-5)
            S_off = QAtlas.fetch(
                TFIM(; J=1.0, h=1.0 + δ), VonNeumannEntropy(), Infinite(); ℓ=ℓ
            )
            # Both finite, same order of magnitude (the gapped CC form
            # tracks the critical log up to a non-universal additive
            # constant that survives the ξ → ∞ limit).
            @test isfinite(S_off)
            @test abs(S_off - S_cr) < 1.0
        end
    end

    @testset "OBC large-N vs Infinite at criticality (sanity check)" begin
        # OBC has a 1/2 prefactor (open boundary halves the universal
        # log coefficient) plus a non-universal boundary entropy.  We
        # only check loose ordering / magnitude, not exact equality.
        model = TFIM(; J=1.0, h=1.0)
        ℓ = 32
        S_obc = QAtlas.fetch(model, VonNeumannEntropy(), OBC(64); ℓ=ℓ)
        S_inf = QAtlas.fetch(model, VonNeumannEntropy(), Infinite(); ℓ=ℓ)
        @test 0 < S_obc
        @test isfinite(S_inf)
        # OBC universal coefficient (c/6) is half the periodic / infinite
        # one (c/3); 2·S_obc ≈ S_inf within an O(1) non-universal offset.
        @test abs(2 * S_obc - S_inf) < 1.5
    end

    @testset "Finite-T at criticality: S(ℓ, β) → S_T=0 as β → ∞" begin
        model = TFIM(; J=1.0, h=1.0)
        ℓ = 50
        S_T0 = QAtlas.fetch(model, VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=Inf)
        S_lowT = QAtlas.fetch(model, VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=1.0e6)
        # (2β/π) sinh(πℓ/β) → 2ℓ as β → ∞ → match T=0 closed form.
        @test S_lowT ≈ S_T0 rtol = 1e-3
    end

    @testset "Finite-T at criticality: high-T linear regime" begin
        # For β ≪ ℓ: sinh(πℓ/β) ≈ exp(πℓ/β)/2 → S ~ (c/3)·(πℓ/β)
        # which is the volume-law thermal entropy of the CFT.
        model = TFIM(; J=1.0, h=1.0)
        c = 0.5
        ℓ = 100
        β = 1.0
        S = QAtlas.fetch(model, VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=β)
        # Asymptote: (c/3) · (πℓ/β + log(β/π))
        S_asymp = (c / 3) * (π * ℓ / β + log(β / π))
        @test S ≈ S_asymp atol = 1e-1
    end

    @testset "Off-critical + finite β: error" begin
        model = TFIM(; J=1.0, h=0.5)
        @test_throws ErrorException QAtlas.fetch(
            model, VonNeumannEntropy(), Infinite(); ℓ=10, beta=1.0
        )
        @test_throws ErrorException QAtlas.fetch(
            model, RenyiEntropy(2.0), Infinite(); ℓ=10, beta=1.0
        )
    end

    @testset "Argument validation: ℓ ≥ 1" begin
        model = TFIM(; J=1.0, h=1.0)
        @test_throws ArgumentError QAtlas.fetch(model, VonNeumannEntropy(), Infinite(); ℓ=0)
        @test_throws ArgumentError QAtlas.fetch(model, RenyiEntropy(2.0), Infinite(); ℓ=0)
    end
end
