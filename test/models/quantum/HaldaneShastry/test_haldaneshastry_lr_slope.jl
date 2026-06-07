using Test
using QAtlas

@testset "HaldaneShastry LR velocity + slope wrappers" begin
    @testset "LiebRobinsonVelocity@Infinite" begin
        m = HaldaneShastry()
        v = QAtlas.fetch(m, QAtlas.LiebRobinsonVelocity(), QAtlas.Infinite())
        @test v ≈ π / 2 atol = 1e-12
        for J in (0.5, 1.0, 2.0, 3.7)
            m2 = HaldaneShastry(; J=J)
            v2 = QAtlas.fetch(m2, QAtlas.LiebRobinsonVelocity(), QAtlas.Infinite())
            @test v2 ≈ π * J / 2 atol = 1e-12
        end
        # negative-J: HaldaneShastry constructor forbids J<=0, test via kwarg override
        @test QAtlas.fetch(
            HaldaneShastry(), QAtlas.LiebRobinsonVelocity(), QAtlas.Infinite(); J=-2.0
        ) ≈ π
    end

    @testset "EntanglementGrowthSlope@Infinite" begin
        m = HaldaneShastry()
        for β in (1.0, 2.0, 5.0, 10.0)
            slope = QAtlas.fetch(
                m, QAtlas.EntanglementGrowthSlope(), QAtlas.Infinite(); beta_eff=β
            )
            @test slope ≈ π^2 / (6 * β) atol = 1e-12
        end
        m2 = HaldaneShastry(; J=2.0)
        slope2 = QAtlas.fetch(
            m2, QAtlas.EntanglementGrowthSlope(), QAtlas.Infinite(); beta_eff=4.0
        )
        @test slope2 ≈ π * 1.0 * π / (3 * 4.0) atol = 1e-12
    end
end
