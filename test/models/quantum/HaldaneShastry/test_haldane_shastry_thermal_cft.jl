# test/models/quantum/HaldaneShastry/test_haldane_shastry_thermal_cft.jl
#
# c=1 CFT low-T stopgap for HaldaneShastry at Infinite() (#524).
# Mirrors test_heisenberg1d_thermal_cft.jl (PR #526) — same v_s, same
# validity gate, different e_0.

using Test
using QAtlas
using QAtlas:
    HaldaneShastry,
    FreeEnergy,
    ThermalEntropy,
    SpecificHeat,
    GroundStateEnergyDensity,
    Infinite

@testset "HaldaneShastry — Infinite c=1 CFT thermal (#524 stopgap)" begin
    @testset "Validity gate at β = 5/J" begin
        m = HaldaneShastry(; J=1.0)
        f_inside = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=5.1)
        f_outside = (@test_logs (:warn,) match_mode=:any QAtlas.fetch(
            m, FreeEnergy(), Infinite(); beta=4.9
        ))
        @test isfinite(f_inside)
        @test isnan(f_outside)
    end

    @testset "LO CFT formula at β = 10/J" begin
        m = HaldaneShastry(; J=1.0)
        f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=10.0)
        e0 = -π^2 / 24
        T = 0.1
        v_s = π / 2
        f_expected = e0 - π * T^2 / (6 * v_s)
        @test isapprox(f, f_expected; atol=1e-12)
    end

    @testset "Entropy and specific heat coincide at LO CFT" begin
        m = HaldaneShastry(; J=1.0)
        s = QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=10.0)
        c = QAtlas.fetch(m, SpecificHeat(), Infinite(); beta=10.0)
        @test isapprox(s, c; atol=1e-12)
        @test isapprox(s, 2 / (3 * 10.0); atol=1e-12)
    end

    @testset "β → ∞ approaches Haldane GS density" begin
        m = HaldaneShastry(; J=1.0)
        e0_inf = QAtlas.fetch(m, GroundStateEnergyDensity(), Infinite())
        f_inf = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1000.0)
        s_inf = QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=1000.0)
        @test isapprox(f_inf, e0_inf; atol=1e-5)
        @test 0 < s_inf < 1e-3
    end

    @testset "f(β) monotone non-decreasing in β" begin
        m = HaldaneShastry(; J=1.0)
        βs = [6.0, 8.0, 12.0, 20.0, 50.0]
        fs = [QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β) for β in βs]
        @test all(diff(fs) .≥ 0)
    end

    @testset "β ≤ 0 raises DomainError" begin
        m = HaldaneShastry(; J=1.0)
        @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=0.0)
        @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=-1.0)
        @test_throws DomainError QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=-1.0)
        @test_throws DomainError QAtlas.fetch(m, SpecificHeat(), Infinite(); beta=-1.0)
    end

    @testset "Scaling f(β, J) / J depends only on β J" begin
        m1 = HaldaneShastry(; J=1.0)
        m2 = HaldaneShastry(; J=2.0)
        f1 = QAtlas.fetch(m1, FreeEnergy(), Infinite(); beta=10.0)
        f2 = QAtlas.fetch(m2, FreeEnergy(), Infinite(); beta=5.0) / 2
        @test isapprox(f1, f2; atol=1e-12)
    end
end
