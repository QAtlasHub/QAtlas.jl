using Test
using QAtlas

@testset "Heisenberg1D entanglement triple" begin
    m = QAtlas.Heisenberg1D()

    @testset "MutualInformation matches Universality(:Heisenberg)" begin
        for (ℓ_A, ℓ_B) in ((5.0, 10.0), (10.0, 20.0), (3.0, 30.0))
            for β in (Inf, 5.0, 1.0)
                mi = QAtlas.fetch(
                    m,
                    QAtlas.MutualInformation(),
                    QAtlas.Infinite();
                    ℓ_A=ℓ_A,
                    ℓ_B=ℓ_B,
                    beta=β,
                )
                ref = QAtlas.fetch(
                    QAtlas.Universality(:Heisenberg),
                    QAtlas.MutualInformation(),
                    QAtlas.Infinite();
                    ℓ_A=ℓ_A,
                    ℓ_B=ℓ_B,
                    beta=β,
                )
                @test mi ≈ ref atol = 1e-12
            end
        end
    end

    @testset "LogarithmicNegativity = (1/4) log[ℓ_A ℓ_B / (ℓ_A + ℓ_B)]" begin
        for (ℓ_A, ℓ_B) in ((5.0, 10.0), (10.0, 20.0))
            e = QAtlas.fetch(
                m, QAtlas.LogarithmicNegativity(), QAtlas.Infinite(); ℓ_A=ℓ_A, ℓ_B=ℓ_B
            )
            ref = (1.0 / 4) * log(ℓ_A * ℓ_B / (ℓ_A + ℓ_B))
            @test e ≈ ref atol = 1e-12
        end
    end

    @testset "EntanglementSaturationDensity = π / (6 β_eff)" begin
        for β in (1.0, 2.0, 5.0, 10.0)
            s = QAtlas.fetch(
                m, QAtlas.EntanglementSaturationDensity(), QAtlas.Infinite(); beta_eff=β
            )
            @test s ≈ π / (6 * β) atol = 1e-12
        end
    end

    @testset "Heisenberg1D ↔ HaldaneShastry agreement (both c = 1)" begin
        m_hs = QAtlas.HaldaneShastry()
        for (ℓ_A, ℓ_B) in ((8.0, 12.0), (5.0, 20.0))
            mi_h1 = QAtlas.fetch(
                m, QAtlas.MutualInformation(), QAtlas.Infinite(); ℓ_A=ℓ_A, ℓ_B=ℓ_B
            )
            mi_hs = QAtlas.fetch(
                m_hs, QAtlas.MutualInformation(), QAtlas.Infinite(); ℓ_A=ℓ_A, ℓ_B=ℓ_B
            )
            @test mi_h1 ≈ mi_hs atol = 1e-12
            ln_h1 = QAtlas.fetch(
                m, QAtlas.LogarithmicNegativity(), QAtlas.Infinite(); ℓ_A=ℓ_A, ℓ_B=ℓ_B
            )
            ln_hs = QAtlas.fetch(
                m_hs, QAtlas.LogarithmicNegativity(), QAtlas.Infinite(); ℓ_A=ℓ_A, ℓ_B=ℓ_B
            )
            @test ln_h1 ≈ ln_hs atol = 1e-12
        end
        for β in (1.0, 5.0)
            s_h1 = QAtlas.fetch(
                m, QAtlas.EntanglementSaturationDensity(), QAtlas.Infinite(); beta_eff=β
            )
            s_hs = QAtlas.fetch(
                m_hs, QAtlas.EntanglementSaturationDensity(), QAtlas.Infinite(); beta_eff=β
            )
            @test s_h1 ≈ s_hs atol = 1e-12
        end
    end
end
