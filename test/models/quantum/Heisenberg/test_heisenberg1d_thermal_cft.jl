# test/models/quantum/Heisenberg/test_heisenberg1d_thermal_cft.jl
#
# c=1 CFT low-T stopgap for Heisenberg1D at Infinite() (#521 Path B).
# Validates: validity gate at β = 5/J, leading-order forms against the
# closed-form values, monotonicity, β → ∞ limit toward Hulthén e₀.

using Test
using QAtlas
using QAtlas: Heisenberg1D, FreeEnergy, ThermalEntropy, SpecificHeat, Infinite

@testset "Heisenberg1D — Infinite c=1 CFT thermal (#521 Path B)" begin

    # ── (a) Validity gate: β ≤ 5/J → NaN + warn ──
    @testset "Validity gate at β = 5/J" begin
        m = Heisenberg1D()
        f_inside = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=5.1, J=1.0)
        f_outside = (@test_logs (:warn,) match_mode=:any QAtlas.fetch(
            m, FreeEnergy(), Infinite(); beta=4.9, J=1.0
        ))
        @test isfinite(f_inside)
        @test isnan(f_outside)
    end

    # ── (b) Leading-order form (β = 10/J) ──
    # f = e₀ - π T² / (6 v_s) = J(1/4 - ln2) - T² / (3J)
    # with T = 0.1, J = 1: f = (0.25 - log 2) - 0.01/3 = -0.4431... - 0.00333...
    @testset "LO CFT formula at β = 10/J" begin
        m = Heisenberg1D()
        J, β = 1.0, 10.0
        f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β, J=J)
        e0 = J * (0.25 - log(2))
        T = 1 / β
        v_s = π * J / 2
        f_expected = e0 - π * T^2 / (6 * v_s)
        @test isapprox(f, f_expected; atol=1e-12)
    end

    @testset "Entropy and specific heat coincide at LO CFT" begin
        m = Heisenberg1D()
        s = QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=10.0, J=1.0)
        c = QAtlas.fetch(m, SpecificHeat(), Infinite(); beta=10.0, J=1.0)
        # At LO CFT s = c = 2T / (3J); check numerical equality.
        @test isapprox(s, c; atol=1e-12)
        @test isapprox(s, 2 / (3 * 10.0); atol=1e-12)
    end

    # ── (c) β → ∞ limit: f → e₀, s → 0, c → 0 ──
    @testset "β → ∞ approaches Hulthén GS density" begin
        m = Heisenberg1D()
        e0 = 0.25 - log(2)  # Hulthén at J=1
        f_inf = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1000.0, J=1.0)
        s_inf = QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=1000.0, J=1.0)
        @test isapprox(f_inf, e0; atol=1e-5)
        @test 0 < s_inf < 1e-3
    end

    # ── (d) Monotonicity: f(β) increases with β at fixed J ──
    @testset "f(β) monotone non-decreasing in β" begin
        m = Heisenberg1D()
        βs = [6.0, 8.0, 12.0, 20.0, 50.0]
        fs = [QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β, J=1.0) for β in βs]
        @test all(diff(fs) .≥ 0)
    end

    # ── (e) β ≤ 0 raises DomainError ──
    @testset "β ≤ 0 raises DomainError" begin
        m = Heisenberg1D()
        @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=0.0, J=1.0)
        @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=-1.0, J=1.0)
        @test_throws DomainError QAtlas.fetch(
            m, ThermalEntropy(), Infinite(); beta=-1.0, J=1.0
        )
        @test_throws DomainError QAtlas.fetch(
            m, SpecificHeat(), Infinite(); beta=-1.0, J=1.0
        )
    end

    # ── (f) J ≤ 0 raises DomainError ──
    @testset "J ≤ 0 raises DomainError" begin
        m = Heisenberg1D()
        @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=10.0, J=0.0)
        @test_throws DomainError QAtlas.fetch(
            m, FreeEnergy(), Infinite(); beta=10.0, J=-1.0
        )
    end

    # ── (g) Scaling: f(β, J) = J · f(β J, 1) by dimensional analysis ──
    @testset "Scaling f(β, J) / J depends only on β J" begin
        m = Heisenberg1D()
        f1 = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=10.0, J=1.0)
        f2 = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=5.0, J=2.0) / 2
        @test isapprox(f1, f2; atol=1e-12)
    end
end
