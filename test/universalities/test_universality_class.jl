using QAtlas
using Test

@testset "UniversalityClass quantity" begin
    # 1. Test fetch for critical quantum model instances
    @testset "Critical quantum models" begin
        # TFIM critical point h = J
        tfim_crit = TFIM(; J=1.0, h=1.0)
        @test fetch(tfim_crit, UniversalityClass(), Infinite()) === Universality(:Ising)

        # XXZ1D critical line (-1 < Δ < 1) and isotropic point (Δ = 1)
        xxz_xy = XXZ1D(; Δ=0.5)
        xxz_heis = XXZ1D(; Δ=1.0)
        @test fetch(xxz_xy, UniversalityClass(), Infinite()) === Universality(:XY)
        @test fetch(xxz_heis, UniversalityClass(), Infinite()) === Universality(:Heisenberg)

        # Heisenberg1D (isotropic AFM point)
        heis1d = Heisenberg1D()
        @test fetch(heis1d, UniversalityClass(), Infinite()) === Universality(:Heisenberg)

        # HaldaneShastry (SU(2)_1 WZW)
        hs = HaldaneShastry()
        @test fetch(hs, UniversalityClass(), Infinite()) === Universality(:Heisenberg)
    end

    # 2. Test fetch for classical model instances at critical point T_c
    @testset "Classical models" begin
        @test fetch(IsingSquare(), UniversalityClass(), Infinite()) === Universality(:Ising)
        @test fetch(IsingTriangular(), UniversalityClass(), Infinite()) === Universality(:Ising)
        @test fetch(CurieWeissIsing(), UniversalityClass(), Infinite()) === Universality(:MeanField)
        @test fetch(TASEP(), UniversalityClass(), Infinite()) === Universality(:KPZ)
    end

    # 3. Test that off-critical model instances throw DomainError
    @testset "Off-critical quantum models throw DomainError" begin
        # TFIM off-critical (h != J)
        tfim_gapped = TFIM(; J=1.0, h=2.0)
        @test_throws DomainError fetch(tfim_gapped, UniversalityClass(), Infinite())

        # XXZ1D off-critical (Δ > 1 or Δ <= -1)
        xxz_gapped = XXZ1D(; Δ=2.0)
        @test_throws DomainError fetch(xxz_gapped, UniversalityClass(), Infinite())
    end

    # 4. Test fetch for Universality{C} objects
    @testset "Universality class objects" begin
        for c in [:Ising, :XY, :Heisenberg, :MeanField, :KPZ]
            u = Universality(c)
            @test fetch(u, UniversalityClass(), Infinite()) === u
        end
    end
end
