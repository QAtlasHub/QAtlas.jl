using Test
using QAtlas

@testset "Universality CFTThermalEntropyDensity (Bloete-Cardy-Nightingale 1986)" begin
    for sym in (:Ising, :Heisenberg, :XY, :Potts3, :Potts4)
        u = QAtlas.Universality(sym)
        c = QAtlas._cardy_central_charge(u)
        for β in (1.0, 2.0, 5.0, 10.0, 0.5)
            s = QAtlas.fetch(
                u, QAtlas.CFTThermalEntropyDensity(), QAtlas.Infinite(); beta=β
            )
            ref = π * Float64(c) / (3 * β)
            @test s ≈ ref atol = 1e-12
        end
    end

    @testset "DomainError on non-positive beta" begin
        u = QAtlas.Universality(:Ising)
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.CFTThermalEntropyDensity(), QAtlas.Infinite(); beta=0.0
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.CFTThermalEntropyDensity(), QAtlas.Infinite(); beta=-2.5
        )
    end

    @testset "Scales as 1/beta" begin
        u = QAtlas.Universality(:Heisenberg)
        s_b1 = QAtlas.fetch(
            u, QAtlas.CFTThermalEntropyDensity(), QAtlas.Infinite(); beta=1.0
        )
        s_b2 = QAtlas.fetch(
            u, QAtlas.CFTThermalEntropyDensity(), QAtlas.Infinite(); beta=2.0
        )
        s_b5 = QAtlas.fetch(
            u, QAtlas.CFTThermalEntropyDensity(), QAtlas.Infinite(); beta=5.0
        )
        @test s_b1 / s_b2 ≈ 2.0 atol = 1e-12
        @test s_b2 / s_b5 ≈ 2.5 atol = 1e-12
    end

    @testset "Thermodynamic relation: e = T s / 2 for CFT" begin
        # e = π c / (6 β²), s = π c / (3 β); T s = s / β = π c / (3 β²) = 2 e
        # so e = T s / 2 ⇔ e / s = T / 2 = 1 / (2β)
        u = QAtlas.Universality(:Ising)
        for β in (1.0, 3.0, 7.0)
            e = QAtlas.fetch(u, QAtlas.ThermalEnergyDensity(), QAtlas.Infinite(); beta=β)
            s = QAtlas.fetch(
                u, QAtlas.CFTThermalEntropyDensity(), QAtlas.Infinite(); beta=β
            )
            @test e / s ≈ 1 / (2 * β) atol = 1e-12
        end
    end
end
