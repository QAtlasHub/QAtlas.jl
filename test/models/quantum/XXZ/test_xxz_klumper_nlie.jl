# test/models/quantum/XXZ/test_xxz_klumper_nlie_dispatch.jl
#
# Part of the Klümper-NLIE-backed FreeEnergy@Infinite test suite for
# `XXZ1D`, split from the original `test_xxz_klumper_nlie.jl` so each
# subset runs in its own CI shard.

using Test
using QAtlas
using QAtlas: XXZ1D, FreeEnergy, ThermalEntropy, SpecificHeat, Energy, Infinite

@testset "XXZ Klümper NLIE — dispatch & validation" begin
    @testset "XX closed-form at Δ = 0" begin
        m = XXZ1D(; J=1.0, Δ=0.0)
        f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1.0)
        @test isfinite(f)
        @test f < 0
    end

    @testset "|Δ| ≥ 1 (gapped) returns NaN + warn" begin
        m = XXZ1D(; J=1.0, Δ=1.5)
        f = (@test_logs (:warn,) match_mode=:any QAtlas.fetch(
            m, FreeEnergy(), Infinite(); beta=1.0
        ))
        @test isnan(f)
    end

    @testset "near-endpoint Δ = -0.999 returns NaN + warn" begin
        m = XXZ1D(; J=1.0, Δ=-0.999)
        f = (@test_logs (:warn,) match_mode=:any QAtlas.fetch(
            m, FreeEnergy(), Infinite(); beta=1.0
        ))
        @test isnan(f)
    end

    @testset "β ≤ 0 raises DomainError" begin
        m = XXZ1D(; J=1.0, Δ=0.5)
        @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=0.0)
        @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=-1.0)
    end

    @testset "J = 0 raises DomainError" begin
        m = XXZ1D(; J=0.0, Δ=0.5)
        @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1.0)
    end

    @testset "build_grid rejects (γ, ε_shift) that would overflow cosh" begin
        γ_unsafe = acos(-0.99)
        @test_throws ArgumentError QAtlas.XXZKlumperNLIE.build_grid(
            γ_unsafe; N=32, L_factor=5.0, ε_shift=0.1
        )
    end

    @testset "ThermalEntropy / SpecificHeat at Δ = 0.5 still NaN (NLIE not yet wired)" begin
        m = XXZ1D(; J=1.0, Δ=0.5)
        s = (@test_logs (:warn,) match_mode=:any QAtlas.fetch(
            m, ThermalEntropy(), Infinite(); beta=1.0
        ))
        c = (@test_logs (:warn,) match_mode=:any QAtlas.fetch(
            m, SpecificHeat(), Infinite(); beta=1.0
        ))
        @test isnan(s)
        @test isnan(c)
    end
end
