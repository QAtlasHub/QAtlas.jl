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

        # Kitaev1D critical line |μ| = 2|t|
        kitaev_crit = Kitaev1D(; μ=2.0, t=1.0, Δ=1.0)
        @test fetch(kitaev_crit, UniversalityClass(), Infinite()) === Universality(:Ising)
    end

    # 2. Test fetch for classical model instances at critical point T_c
    @testset "Classical models" begin
        @test fetch(IsingSquare(), UniversalityClass(), Infinite()) === Universality(:Ising)
        @test fetch(IsingTriangular(), UniversalityClass(), Infinite()) === Universality(:Ising)
        @test fetch(CurieWeissIsing(), UniversalityClass(), Infinite()) === Universality(:MeanField)
        @test fetch(TASEP(), UniversalityClass(), Infinite()) === Universality(:KPZ)

        # ZnClock
        @test fetch(ZnClock(; n=2), UniversalityClass(), Infinite()) === Universality(:Ising)
        @test fetch(ZnClock(; n=3), UniversalityClass(), Infinite()) === Universality(:Potts3)

        # ZnParafermion
        @test fetch(ZnParafermion(; n=2), UniversalityClass(), Infinite()) === Universality(:Ising)
        @test fetch(ZnParafermion(; n=3), UniversalityClass(), Infinite()) === Universality(:Potts3)
        @test fetch(ZnParafermion(; n=4), UniversalityClass(), Infinite()) === Universality(:Potts4)

        # SixVertex (disordered phase |Δ| < 1)
        @test fetch(SixVertex(; a=1.0, b=1.0, c=1.0), UniversalityClass(), Infinite()) === Universality(:XY)

        # DimerLattice (c = 1 compact free boson)
        @test fetch(DimerLattice(), UniversalityClass(), Infinite()) === Universality(:XY)

        # TricriticalIsing
        @test fetch(TricriticalIsing(), UniversalityClass(), Infinite()) === Universality(:TricriticalIsing)

        # TricriticalPotts3
        @test fetch(TricriticalPotts3(), UniversalityClass(), Infinite()) === Universality(:TricriticalPotts3)
    end

    # 3. Test that off-critical model instances throw DomainError
    @testset "Off-critical quantum models throw DomainError" begin
        # TFIM off-critical (h != J)
        tfim_gapped = TFIM(; J=1.0, h=2.0)
        @test_throws DomainError fetch(tfim_gapped, UniversalityClass(), Infinite())

        # XXZ1D off-critical (Δ > 1 or Δ <= -1)
        xxz_gapped = XXZ1D(; Δ=2.0)
        @test_throws DomainError fetch(xxz_gapped, UniversalityClass(), Infinite())

        # Kitaev1D off-critical (|μ| != 2|t|)
        kitaev_gapped = Kitaev1D(; μ=3.0, t=1.0, Δ=1.0)
        @test_throws DomainError fetch(kitaev_gapped, UniversalityClass(), Infinite())

        # ZnClock unsupported / off-critical (n >= 4)
        @test_throws DomainError fetch(ZnClock(; n=4), UniversalityClass(), Infinite())

        # ZnParafermion unsupported / off-critical (n = 5)
        @test_throws DomainError fetch(ZnParafermion(; n=5), UniversalityClass(), Infinite())

        # SixVertex off-critical (|Δ| >= 1, e.g. ferroelectric or AFE)
        @test_throws DomainError fetch(SixVertex(; a=3.0, b=1.0, c=1.0), UniversalityClass(), Infinite())
    end

    # 4. Test fetch for Universality{C} objects
    @testset "Universality class objects" begin
        for c in [:Ising, :XY, :Heisenberg, :MeanField, :KPZ, :Potts3, :Potts4, :TricriticalIsing, :TricriticalPotts3]
            u = Universality(c)
            @test fetch(u, UniversalityClass(), Infinite()) === u
        end
    end
end

