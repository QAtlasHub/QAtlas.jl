using Test
using QAtlas

@testset "XXZ1D LR velocity + slope at XX point (Δ = 0)" begin
    @testset "LiebRobinsonVelocity = 2|J|" begin
        for J in (0.5, 1.0, 2.0, 5.0)
            m = QAtlas.XXZ1D(; J=J, Δ=0.0)
            v = QAtlas.fetch(m, QAtlas.LiebRobinsonVelocity(), QAtlas.Infinite())
            @test v ≈ 2 * abs(J) atol = 1e-12
        end
        m_neg = QAtlas.XXZ1D(; J=1.0, Δ=0.0)
        v_neg = QAtlas.fetch(
            m_neg, QAtlas.LiebRobinsonVelocity(), QAtlas.Infinite(); J=-3.7
        )
        @test v_neg ≈ 2 * 3.7 atol = 1e-12
    end

    @testset "EntanglementGrowthSlope = πcv/(3β_eff) with c=1, v=2|J|" begin
        for J in (1.0, 2.0)
            m = QAtlas.XXZ1D(; J=J, Δ=0.0)
            for β in (1.0, 2.0, 5.0, 10.0)
                slope = QAtlas.fetch(
                    m, QAtlas.EntanglementGrowthSlope(), QAtlas.Infinite(); beta_eff=β
                )
                @test slope ≈ 2π * abs(J) / (3 * β) atol = 1e-12
            end
        end
    end

    @testset "DomainError off the XX point (Δ != 0)" begin
        for Δ in (0.1, 0.5, 1.0, -0.5)
            m = QAtlas.XXZ1D(; J=1.0, Δ=Δ)
            @test_throws DomainError QAtlas.fetch(
                m, QAtlas.LiebRobinsonVelocity(), QAtlas.Infinite()
            )
            @test_throws DomainError QAtlas.fetch(
                m, QAtlas.EntanglementGrowthSlope(), QAtlas.Infinite(); beta_eff=1.0
            )
        end
    end
end
