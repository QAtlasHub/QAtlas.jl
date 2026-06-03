using Test
using QAtlas
using QAtlas: Universality, EntanglementGrowthSlope, Infinite, CentralCharge

@testset "Universality EntanglementGrowthSlope (Calabrese-Cardy 2005, #580)" begin
    @testset "Formula πcv/(3β_eff) across 5 universality classes" begin
        d_for = Dict(:Ising=>2, :XY=>2, :Heisenberg=>1, :Potts3=>2, :Potts4=>2)
        for class in (:Ising, :XY, :Heisenberg, :Potts3, :Potts4)
            c = float(QAtlas.fetch(Universality(class), CentralCharge(); d=d_for[class]))
            for v in (1.0, 2.0, 5.0), β_eff in (1.0, 5.0, 50.0)
                slope = QAtlas.fetch(
                    Universality(class),
                    EntanglementGrowthSlope(),
                    Infinite();
                    v=v,
                    beta_eff=β_eff,
                )
                expected = π * c * v / (3 * β_eff)
                @test slope ≈ expected atol = 1e-14
            end
        end
    end

    @testset "Linearity in v (fixed c, β_eff)" begin
        for class in (:Ising, :Heisenberg)
            d = class == :Heisenberg ? 1 : 2
            β_eff = 4.0
            s1 = QAtlas.fetch(
                Universality(class),
                EntanglementGrowthSlope(),
                Infinite();
                v=1.0,
                beta_eff=β_eff,
            )
            s3 = QAtlas.fetch(
                Universality(class),
                EntanglementGrowthSlope(),
                Infinite();
                v=3.0,
                beta_eff=β_eff,
            )
            @test s3 ≈ 3 * s1 atol = 1e-14
        end
    end

    @testset "Inverse linearity in β_eff (fixed c, v)" begin
        v = 2.5
        for class in (:Ising, :Heisenberg)
            s1 = QAtlas.fetch(
                Universality(class),
                EntanglementGrowthSlope(),
                Infinite();
                v=v,
                beta_eff=1.0,
            )
            s10 = QAtlas.fetch(
                Universality(class),
                EntanglementGrowthSlope(),
                Infinite();
                v=v,
                beta_eff=10.0,
            )
            @test s10 ≈ s1 / 10 atol = 1e-14
        end
    end

    @testset "Argument validation: v > 0, β_eff > 0" begin
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising),
            EntanglementGrowthSlope(),
            Infinite();
            v=-1.0,
            beta_eff=1.0,
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), EntanglementGrowthSlope(), Infinite(); v=1.0, beta_eff=0.0
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Heisenberg),
            EntanglementGrowthSlope(),
            Infinite();
            v=0.0,
            beta_eff=1.0,
        )
    end
end
