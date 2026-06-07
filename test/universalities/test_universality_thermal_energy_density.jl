using Test
using QAtlas

@testset "Universality ThermalEnergyDensity (Affleck 1986)" begin
    for sym in (:Ising, :Heisenberg, :XY, :Potts3, :Potts4)
        u = QAtlas.Universality(sym)
        c = QAtlas._cardy_central_charge(u)
        for β in (1.0, 2.0, 5.0, 10.0, 0.5)
            e_th = QAtlas.fetch(u, QAtlas.ThermalEnergyDensity(), QAtlas.Infinite(); beta=β)
            ref = π * Float64(c) / (6 * β^2)
            @test e_th ≈ ref atol = 1e-12
        end
    end

    @testset "DomainError on non-positive beta" begin
        u = QAtlas.Universality(:Ising)
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.ThermalEnergyDensity(), QAtlas.Infinite(); beta=0.0
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.ThermalEnergyDensity(), QAtlas.Infinite(); beta=-1.0
        )
    end

    @testset "Linear in c (Ising vs Heisenberg ratio = 1/2)" begin
        β = 3.0
        u_is = QAtlas.Universality(:Ising)
        u_he = QAtlas.Universality(:Heisenberg)
        e_is = QAtlas.fetch(u_is, QAtlas.ThermalEnergyDensity(), QAtlas.Infinite(); beta=β)
        e_he = QAtlas.fetch(u_he, QAtlas.ThermalEnergyDensity(), QAtlas.Infinite(); beta=β)
        @test e_is / e_he ≈ 0.5 atol = 1e-12
    end

    @testset "Scales as 1/beta^2" begin
        u = QAtlas.Universality(:Heisenberg)
        e_b1 = QAtlas.fetch(u, QAtlas.ThermalEnergyDensity(), QAtlas.Infinite(); beta=1.0)
        e_b2 = QAtlas.fetch(u, QAtlas.ThermalEnergyDensity(), QAtlas.Infinite(); beta=2.0)
        e_b4 = QAtlas.fetch(u, QAtlas.ThermalEnergyDensity(), QAtlas.Infinite(); beta=4.0)
        @test e_b1 / e_b2 ≈ 4.0 atol = 1e-12
        @test e_b2 / e_b4 ≈ 4.0 atol = 1e-12
    end
end
