using Test
using QAtlas
using QAtlas: Heisenberg1D, VonNeumannEntropy, RenyiEntropy, Infinite

@testset "Heisenberg1D Calabrese-Cardy entanglement at Infinite (#580 Phase 1)" begin
    c = 1.0  # SU(2) Heisenberg point sits on c=1 free-boson CFT

    @testset "T = 0 critical: VN = (c/3) log(2 ℓ)" begin
        m = Heisenberg1D()
        for ℓ in (1, 5, 10, 50, 100)
            S = QAtlas.fetch(m, VonNeumannEntropy(), Infinite(); ℓ=ℓ)
            @test S ≈ (c / 3) * log(2 * ℓ) atol = 1e-12
        end
    end

    @testset "T = 0 critical: Renyi-α coefficient (c/6)(1 + 1/α) log(2 ℓ)" begin
        m = Heisenberg1D()
        for α in (0.5, 2.0, 3.0, 5.0), ℓ in (5, 20, 50)
            S = QAtlas.fetch(m, RenyiEntropy(α), Infinite(); ℓ=ℓ)
            @test S ≈ (c / 6) * (1 + 1 / α) * log(2 * ℓ) atol = 1e-12
        end
    end

    @testset "Finite T critical: VN matches CC sinh form" begin
        m = Heisenberg1D()
        for ℓ in (5, 20), β in (1.0, 5.0, 20.0)
            S = QAtlas.fetch(m, VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=β)
            expected = (c / 3) * log((2 * β / π) * sinh(π * ℓ / β))
            @test S ≈ expected atol = 1e-12
        end
    end

    @testset "β → ∞ recovery: S(ℓ, β=large) ≈ S_T=0" begin
        m = Heisenberg1D()
        ℓ = 10
        S_T0 = QAtlas.fetch(m, VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=Inf)
        S_lowT = QAtlas.fetch(m, VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=1.0e6)
        # (2β/π) sinh(πℓ/β) → 2ℓ as β → ∞ ⇒ matches T=0 closed form.
        @test S_lowT ≈ S_T0 rtol = 1e-3
    end

    @testset "High-T linear regime: S ~ (c/3) · (πℓ/β + log(β/π))" begin
        # For β ≪ ℓ: sinh(πℓ/β) ≈ exp(πℓ/β)/2 ⇒ S ~ (c/3) · (πℓ/β + log(β/π)).
        m = Heisenberg1D()
        ℓ = 20
        β = 1.0
        S = QAtlas.fetch(m, VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=β)
        S_asymp = (c / 3) * (π * ℓ / β + log(β / π))
        @test S ≈ S_asymp atol = 1e-1
    end

    @testset "Argument errors" begin
        m = Heisenberg1D()
        @test_throws ArgumentError QAtlas.fetch(m, VonNeumannEntropy(), Infinite(); ℓ=0)
        @test_throws ArgumentError QAtlas.fetch(m, RenyiEntropy(2.0), Infinite(); ℓ=-1)
    end
end
