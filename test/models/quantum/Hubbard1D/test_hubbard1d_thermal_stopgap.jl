# test/models/quantum/Hubbard1D/test_hubbard1d_thermal_stopgap.jl
#
# Phase-2A regime-based delegation tests (#523 stopgap).
# Covers the 4 regimes:
#   (a) U = 0 exact -> 2 x TightBinding1D
#   (b) very-high T (β·max(t,U) <= 0.05) -> -T ln 4
#   (c) strong-coupling + low T (U/t>=10, β·J_eff>=5) -> Lieb-Wu e_0 + Heisenberg CFT excess
#   (d) intermediate -> NaN + warn
# Plus validation: half-filling guard, beta/t/U DomainError.

using Test
using QAtlas
using QAtlas:
    Hubbard1D,
    TightBinding1D,
    FreeEnergy,
    ThermalEntropy,
    SpecificHeat,
    GroundStateEnergyDensity,
    Infinite

@testset "Hubbard1D — finite-T Phase-2A stopgap (#523)" begin
    @testset "(A) U = 0 delegates to 2 x TightBinding1D" begin
        m = Hubbard1D(1.0, 0.0, 0.0)  # t=1, U=0, μ=U/2=0
        tb = TightBinding1D(; t=1.0, μ=0.0)
        β = 1.0
        f_h = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
        f_tb = QAtlas.fetch(tb, FreeEnergy(), Infinite(); beta=β)
        @test isapprox(f_h, 2 * f_tb; atol=1e-10)

        s_h = QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=β)
        s_tb = QAtlas.fetch(tb, ThermalEntropy(), Infinite(); beta=β)
        @test isapprox(s_h, 2 * s_tb; atol=1e-10)
    end

    @testset "(B) very-high-T: -T ln 4 / ln 4 / 0" begin
        m = Hubbard1D(1.0, 4.0, 2.0)  # half-filling
        β = 0.01  # β·max(t,U) = 0.04 < 0.05
        f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
        @test isapprox(f, -log(4) / β; atol=1e-10)
        @test isapprox(
            QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=β), log(4); atol=1e-10
        )
        @test isapprox(QAtlas.fetch(m, SpecificHeat(), Infinite(); beta=β), 0.0; atol=1e-10)
    end

    @testset "(C) strong-coupling + low-T: e_0 - T^2/(3 J_eff)" begin
        t, U = 1.0, 20.0  # U/t = 20 >= 10
        m = Hubbard1D(t, U, U/2)
        J_eff = 4 * t^2 / U  # = 0.2
        β = 50.0  # β·J_eff = 10 >= 5
        f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
        e0 = QAtlas.fetch(m, GroundStateEnergyDensity(), Infinite())
        T = 1/β
        @test isapprox(f, e0 - T^2 / (3 * J_eff); atol=1e-12)

        # Entropy & specific heat in this regime
        s = QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=β)
        c = QAtlas.fetch(m, SpecificHeat(), Infinite(); beta=β)
        @test isapprox(s, 2 * T / (3 * J_eff); atol=1e-12)
        @test isapprox(c, 2 * T / (3 * J_eff); atol=1e-12)
    end

    @testset "(D) intermediate regime: NaN + warn" begin
        m = Hubbard1D(1.0, 4.0, 2.0)  # U/t = 4 < 10
        β = 1.0  # β·max(t,U) = 4 >> 0.05
        f = (@test_logs (:warn,) match_mode=:any QAtlas.fetch(
            m, FreeEnergy(), Infinite(); beta=β
        ))
        @test isnan(f)

        s = (@test_logs (:warn,) match_mode=:any QAtlas.fetch(
            m, ThermalEntropy(), Infinite(); beta=β
        ))
        c = (@test_logs (:warn,) match_mode=:any QAtlas.fetch(
            m, SpecificHeat(), Infinite(); beta=β
        ))
        @test isnan(s)
        @test isnan(c)
    end

    @testset "Off half-filling raises DomainError" begin
        m = Hubbard1D(1.0, 4.0, 1.0)  # μ = 1, but U/2 = 2 -> off half-filling
        @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1.0)
        @test_throws DomainError QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=1.0)
        @test_throws DomainError QAtlas.fetch(m, SpecificHeat(), Infinite(); beta=1.0)
    end

    @testset "β <= 0 raises DomainError" begin
        m = Hubbard1D(1.0, 4.0, 2.0)
        @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=0.0)
        @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=-1.0)
    end

    @testset "C regime f(β) monotone non-decreasing in β" begin
        t, U = 1.0, 20.0
        m = Hubbard1D(t, U, U/2)
        # β·J_eff: 50->10, 100->20, both >= 5
        β_ok = [50.0, 100.0]
        fs = [QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β) for β in β_ok]
        @test all(diff(fs) .≥ 0)
    end
end
