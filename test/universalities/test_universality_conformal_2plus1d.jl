# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: 2+1D CFT universalities (SphereFreeEnergy & CornerEntanglementCoefficient)
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Universality :: 2+1D CFT universalities" begin
    # ── 1. SphereFreeEnergy ──
    @testset "SphereFreeEnergy values in d=3" begin
        # Ising
        @test QAtlas.fetch(Universality(:Ising), SphereFreeEnergy(), Infinite(); d=3) ≈ 0.0612 atol=1e-6
        # XY
        @test QAtlas.fetch(Universality(:XY), SphereFreeEnergy(), Infinite(); d=3) ≈ 0.121 atol=1e-6
        # Heisenberg
        @test QAtlas.fetch(Universality(:Heisenberg), SphereFreeEnergy(), Infinite(); d=3) ≈ 0.180 atol=1e-6
    end

    @testset "SphereFreeEnergy dimension guard" begin
        @test_throws ArgumentError QAtlas.fetch(Universality(:Ising), SphereFreeEnergy(), Infinite(); d=2)
        @test_throws ArgumentError QAtlas.fetch(Universality(:Ising), SphereFreeEnergy(), Infinite(); d=4)
    end

    @testset "SphereFreeEnergy unsupported classes" begin
        @test_throws ArgumentError QAtlas.fetch(Universality(:Potts3), SphereFreeEnergy(), Infinite(); d=3)
    end

    # ── 2. CornerEntanglementCoefficient ──
    @testset "CornerEntanglementCoefficient prefactor σ" begin
        # Ising: σ ≈ 0.0036974
        @test QAtlas.fetch(Universality(:Ising), CornerEntanglementCoefficient(), Infinite(); d=3) ≈ 0.0036974 atol=1e-7
        # XY: σ ≈ 0.0073773
        @test QAtlas.fetch(Universality(:XY), CornerEntanglementCoefficient(), Infinite(); d=3) ≈ 0.0073773 atol=1e-7
        # Heisenberg: σ ≈ 0.011050
        @test QAtlas.fetch(Universality(:Heisenberg), CornerEntanglementCoefficient(), Infinite(); d=3) ≈ 0.011050 atol=1e-7
    end

    @testset "CornerEntanglementCoefficient angle-dependent c(θ)" begin
        theta = 2.0  # 0 ≤ theta ≤ π
        # c(θ) ≈ σ * (π - θ)^2
        # Ising
        σ_ising = 0.0036974
        @test QAtlas.fetch(Universality(:Ising), CornerEntanglementCoefficient(), Infinite(); d=3, theta=theta) ≈ σ_ising * (π - theta)^2 atol=1e-10

        # XY
        σ_xy = 0.0073773
        @test QAtlas.fetch(Universality(:XY), CornerEntanglementCoefficient(), Infinite(); d=3, theta=theta) ≈ σ_xy * (π - theta)^2 atol=1e-10
    end

    @testset "CornerEntanglementCoefficient validations" begin
        # Dimension guard
        @test_throws ArgumentError QAtlas.fetch(Universality(:Ising), CornerEntanglementCoefficient(), Infinite(); d=2)
        # Theta boundary guards
        @test_throws ArgumentError QAtlas.fetch(Universality(:Ising), CornerEntanglementCoefficient(), Infinite(); d=3, theta=-0.1)
        @test_throws ArgumentError QAtlas.fetch(Universality(:Ising), CornerEntanglementCoefficient(), Infinite(); d=3, theta=4.0)
    end
end
