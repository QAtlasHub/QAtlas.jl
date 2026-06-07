using Test
using QAtlas

@testset "Heisenberg1D LR velocity + slope (des Cloizeaux-Pearson)" begin
    m = QAtlas.Heisenberg1D()

    @testset "LiebRobinsonVelocity = π J / 2" begin
        for J in (0.5, 1.0, 2.0, 5.0)
            v = QAtlas.fetch(m, QAtlas.LiebRobinsonVelocity(), QAtlas.Infinite(); J=J)
            @test v ≈ π * J / 2 atol = 1e-12
        end
        v_neg = QAtlas.fetch(m, QAtlas.LiebRobinsonVelocity(), QAtlas.Infinite(); J=-3.0)
        @test v_neg ≈ 3π / 2 atol = 1e-12
        v_default = QAtlas.fetch(m, QAtlas.LiebRobinsonVelocity(), QAtlas.Infinite())
        @test v_default ≈ π / 2 atol = 1e-12
    end

    @testset "EntanglementGrowthSlope = π² J / (6 β_eff)" begin
        for β in (1.0, 2.0, 5.0)
            slope = QAtlas.fetch(
                m, QAtlas.EntanglementGrowthSlope(), QAtlas.Infinite(); beta_eff=β
            )
            @test slope ≈ π^2 / (6 * β) atol = 1e-12
        end
    end

    @testset "Heisenberg1D ↔ HaldaneShastry LR velocity equivalence at default J" begin
        m_hs = QAtlas.HaldaneShastry()
        v_h1 = QAtlas.fetch(m, QAtlas.LiebRobinsonVelocity(), QAtlas.Infinite())
        v_hs = QAtlas.fetch(m_hs, QAtlas.LiebRobinsonVelocity(), QAtlas.Infinite())
        @test v_h1 ≈ v_hs atol = 1e-12
    end
end
