using Test
using QAtlas

@testset "Three-model EntanglementSaturationDensity wrappers" begin
    @testset "TFIM h = J critical → Universality(:Ising) (c = 1/2)" begin
        for J in (1.0, 2.0)
            m = QAtlas.TFIM(; h=J, J=J)
            for β in (1.0, 2.0, 5.0)
                s = QAtlas.fetch(
                    m, QAtlas.EntanglementSaturationDensity(), QAtlas.Infinite(); beta_eff=β
                )
                ref = π * 0.5 / (6 * β)
                @test s ≈ ref atol = 1e-12
            end
        end
        m_off = QAtlas.TFIM(; h=0.5, J=1.0)
        @test_throws DomainError QAtlas.fetch(
            m_off, QAtlas.EntanglementSaturationDensity(), QAtlas.Infinite(); beta_eff=1.0
        )
    end

    @testset "HaldaneShastry → Universality(:Heisenberg) (c = 1)" begin
        for J in (0.5, 1.0, 2.0)
            m = QAtlas.HaldaneShastry(; J=J)
            for β in (1.0, 2.0, 5.0)
                s = QAtlas.fetch(
                    m, QAtlas.EntanglementSaturationDensity(), QAtlas.Infinite(); beta_eff=β
                )
                ref = π / (6 * β)
                @test s ≈ ref atol = 1e-12
            end
        end
    end

    @testset "XXZ1D at XX (Δ = 0) → Universality(:XY) (c = 1)" begin
        for J in (0.5, 1.0, 2.0)
            m = QAtlas.XXZ1D(; J=J, Δ=0.0)
            for β in (1.0, 2.0, 5.0)
                s = QAtlas.fetch(
                    m, QAtlas.EntanglementSaturationDensity(), QAtlas.Infinite(); beta_eff=β
                )
                ref = π / (6 * β)
                @test s ≈ ref atol = 1e-12
            end
        end
        for Δ in (0.5, 1.0, -0.5)
            m_off = QAtlas.XXZ1D(; J=1.0, Δ=Δ)
            @test_throws DomainError QAtlas.fetch(
                m_off,
                QAtlas.EntanglementSaturationDensity(),
                QAtlas.Infinite();
                beta_eff=1.0,
            )
        end
    end

    @testset "Universality central-charge ratio reflected in 3 models" begin
        β = 3.0
        m_tfim = QAtlas.TFIM(; h=1.0, J=1.0)
        m_hs = QAtlas.HaldaneShastry(; J=1.0)
        s_tfim = QAtlas.fetch(
            m_tfim, QAtlas.EntanglementSaturationDensity(), QAtlas.Infinite(); beta_eff=β
        )
        s_hs = QAtlas.fetch(
            m_hs, QAtlas.EntanglementSaturationDensity(), QAtlas.Infinite(); beta_eff=β
        )
        @test s_tfim / s_hs ≈ 0.5 atol = 1e-12
    end
end
