# test/models/quantum/HaldaneShastry/test_haldane_shastry.jl
#
# Scaffold tests for the spin-1/2 Haldane-Shastry chain (#225):
#   (a) GroundStateEnergyDensity@Infinite closed form E_0/N = -π² J / 24
#   (b) spinon dispersion ε(0) = ε(π) = 0, max at k = π/2 equals J π² / 8
#   (c) sound velocity v_s = π J / 2 (matches Heisenberg c=1 SU(2)_1)
#   (d) construction validation: J ≤ 0 raises DomainError
#   (e) dispersion outside [0, π] raises DomainError
#   (f) J scaling: E_0/N linear in J

using Test
using QAtlas
using QAtlas:
    HaldaneShastry,
    GroundStateEnergyDensity,
    Infinite,
    haldane_shastry_spinon_dispersion,
    haldane_shastry_sound_velocity

@testset "HaldaneShastry — scaffold (#225)" begin
    @testset "GS energy density E_0/N = -π² J / 24" begin
        m = HaldaneShastry()
        e0 = QAtlas.fetch(m, GroundStateEnergyDensity(), Infinite())
        @test isapprox(e0, -π^2 / 24; atol=1e-12)
    end

    @testset "Spinon dispersion endpoints + maximum" begin
        m = HaldaneShastry(; J=1.0)
        @test haldane_shastry_spinon_dispersion(m, 0.0) ≈ 0 atol=1e-12
        @test haldane_shastry_spinon_dispersion(m, Float64(π)) ≈ 0 atol=1e-12
        @test isapprox(haldane_shastry_spinon_dispersion(m, π / 2), π^2 / 8; atol=1e-12)
    end

    @testset "Sound velocity v_s = π J / 2 (Heisenberg-class CFT)" begin
        m = HaldaneShastry(; J=2.0)
        @test isapprox(haldane_shastry_sound_velocity(m), π; atol=1e-12)
    end

    @testset "J ≤ 0 raises DomainError" begin
        @test_throws DomainError HaldaneShastry(; J=0.0)
        @test_throws DomainError HaldaneShastry(; J=-1.0)
    end

    @testset "Dispersion outside [0, π] raises DomainError" begin
        m = HaldaneShastry()
        @test_throws DomainError haldane_shastry_spinon_dispersion(m, -0.1)
        @test_throws DomainError haldane_shastry_spinon_dispersion(m, π + 0.1)
    end

    @testset "E_0/N linear in J" begin
        e0_1 = QAtlas.fetch(HaldaneShastry(; J=1.0), GroundStateEnergyDensity(), Infinite())
        e0_3 = QAtlas.fetch(HaldaneShastry(; J=3.0), GroundStateEnergyDensity(), Infinite())
        @test isapprox(e0_3, 3 * e0_1; atol=1e-12)
    end
end
