using Test
using QAtlas
using QAtlas: TFIM, EntanglementGrowthSlope, Infinite, Universality, LiebRobinsonVelocity

@testset "TFIM EntanglementGrowthSlope wrapper at h=J critical (#580 + #579 cross)" begin
    @testset "Slope = π J / (3 β_eff) at h = J critical" begin
        for J in (0.5, 1.0, 2.0, 5.0), β in (1.0, 5.0, 20.0)
            m = TFIM(; J=J, h=J)
            slope = QAtlas.fetch(m, EntanglementGrowthSlope(), Infinite(); beta_eff=β)
            expected = π * J / (3 * β)
            @test slope ≈ expected atol = 1e-14
        end
    end

    @testset "End-to-end chain: c(Ising) · v(TFIM) / (3 β_eff)" begin
        m = TFIM(; J=1.5, h=1.5)
        β = 4.0
        slope_model = QAtlas.fetch(m, EntanglementGrowthSlope(), Infinite(); beta_eff=β)
        c = float(QAtlas.fetch(Universality(:Ising), QAtlas.CentralCharge(); d=2))
        v = QAtlas.fetch(m, LiebRobinsonVelocity(), Infinite())
        slope_chain = π * c * v / (3 * β)
        @test slope_model ≈ slope_chain atol = 1e-14
        # Verify the values: c = 1/2, v = 2J = 3.0, expected = π·0.5·3 / 12 = π/8
        @test slope_chain ≈ π / 8 atol = 1e-14
    end

    @testset "Off-critical h ≠ J: DomainError" begin
        for (J, h) in ((1.0, 0.5), (1.0, 2.0), (1.0, 1.5))
            m = TFIM(; J=J, h=h)
            @test_throws DomainError QAtlas.fetch(
                m, EntanglementGrowthSlope(), Infinite(); beta_eff=5.0
            )
        end
    end

    @testset "Critical detection tolerance" begin
        # Within 1e-10 of criticality should still delegate (not throw).
        m = TFIM(; J=1.0, h=1.0 + 1e-12)
        slope = QAtlas.fetch(m, EntanglementGrowthSlope(), Infinite(); beta_eff=1.0)
        @test isfinite(slope)
        @test slope > 0
    end
end
