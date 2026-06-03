using Test
using QAtlas
using QAtlas:
    Universality,
    EntanglementSaturationDensity,
    EntanglementGrowthSlope,
    Infinite,
    CentralCharge

@testset "Universality EntanglementSaturationDensity = π c / (6 β_eff) CC 2005 (#580)" begin
    @testset "Closed form across 5 universality classes" begin
        d_for = Dict(:Ising=>2, :XY=>2, :Heisenberg=>1, :Potts3=>2, :Potts4=>2)
        for class in (:Ising, :XY, :Heisenberg, :Potts3, :Potts4)
            c = float(QAtlas.fetch(Universality(class), CentralCharge(); d=d_for[class]))
            for β in (0.5, 1.0, 5.0, 20.0, 100.0)
                s = QAtlas.fetch(
                    Universality(class),
                    EntanglementSaturationDensity(),
                    Infinite();
                    beta_eff=β,
                )
                expected = π * c / (6 * β)
                @test s ≈ expected atol = 1e-14
            end
        end
    end

    @testset "Crossover relation: s_sat * (2 v) = slope * 1 (per universality + β_eff)" begin
        # CC linear-regime extrapolation matches saturation at t = L / (2 v):
        # slope * (L / (2 v)) = saturation * L  =>  s_sat / slope = 1 / (2 v)
        for class in (:Ising, :Heisenberg)
            for v in (1.0, 2.0, 5.0), β in (1.0, 10.0)
                s_sat = QAtlas.fetch(
                    Universality(class),
                    EntanglementSaturationDensity(),
                    Infinite();
                    beta_eff=β,
                )
                slope = QAtlas.fetch(
                    Universality(class),
                    EntanglementGrowthSlope(),
                    Infinite();
                    v=v,
                    beta_eff=β,
                )
                @test s_sat / slope ≈ 1 / (2 * v) atol = 1e-14
            end
        end
    end

    @testset "Inverse linearity in β_eff at fixed c" begin
        for class in (:Ising, :Heisenberg)
            s1 = QAtlas.fetch(
                Universality(class),
                EntanglementSaturationDensity(),
                Infinite();
                beta_eff=1.0,
            )
            s10 = QAtlas.fetch(
                Universality(class),
                EntanglementSaturationDensity(),
                Infinite();
                beta_eff=10.0,
            )
            @test s10 ≈ s1 / 10 atol = 1e-14
        end
    end

    @testset "Linearity in c at fixed β_eff" begin
        β = 4.0
        s_Ising = QAtlas.fetch(
            Universality(:Ising), EntanglementSaturationDensity(), Infinite(); beta_eff=β
        )
        s_Heis = QAtlas.fetch(
            Universality(:Heisenberg),
            EntanglementSaturationDensity(),
            Infinite();
            beta_eff=β,
        )
        c_I = float(QAtlas.fetch(Universality(:Ising), CentralCharge(); d=2))
        c_H = float(QAtlas.fetch(Universality(:Heisenberg), CentralCharge(); d=1))
        @test s_Heis / s_Ising ≈ c_H / c_I atol = 1e-12
    end

    @testset "Argument validation" begin
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), EntanglementSaturationDensity(), Infinite(); beta_eff=-1.0
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Heisenberg),
            EntanglementSaturationDensity(),
            Infinite();
            beta_eff=0.0,
        )
    end
end
