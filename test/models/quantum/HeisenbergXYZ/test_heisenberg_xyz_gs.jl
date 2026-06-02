# test/models/quantum/HeisenbergXYZ/test_heisenberg_xyz_gs.jl
#
# Phase-2 ground-state energy density coverage for HeisenbergXYZ
# (XY anisotropic line via Lieb-Schultz-Mattis 1961 + axial XXZ delegation).

using Test
using QuadGK: quadgk
using QAtlas
using QAtlas: fetch, GroundStateEnergyDensity, Infinite, HeisenbergXYZ

@testset "HeisenbergXYZ — GroundStateEnergyDensity@Infinite (Phase 2)" begin
    @testset "Heisenberg isotropic limit (Jx=Jy=Jz=J, AFM)" begin
        # Bethe 1931 / Hulthen 1938: ε₀ = J·(1/4 - ln 2)
        m = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=1.0)
        f = fetch(m, GroundStateEnergyDensity(), Infinite())
        @test f ≈ 0.25 - log(2) atol=1e-12
    end

    @testset "XX free-fermion point (Jx=Jy=1, Jz=0) — XXZ delegation" begin
        # ε₀ = -1/π (free-fermion XX chain, QAtlas S=σ/2 normalization)
        m = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=0.0)
        f = fetch(m, GroundStateEnergyDensity(), Infinite())
        @test f ≈ -1/π atol=1e-12
    end

    @testset "Anisotropic XY line (Jx=2, Jy=1, Jz=0) — LSM closed form" begin
        # Strict anisotropy on Jz=0 line: LSM 1961 integral
        m = HeisenbergXYZ(; Jx=2.0, Jy=1.0, Jz=0.0)
        f = fetch(m, GroundStateEnergyDensity(), Infinite())
        # Verify against direct LSM integral evaluation
        ref, _ = quadgk(k -> sqrt(4 + 1 + 4 * cos(2k)), 0.0, π / 2; rtol=1e-12)
        @test f ≈ -ref / (4π) atol=1e-10
        # Sanity bound: should lie between Heisenberg (-0.443) and 0
        @test -0.5 < f < -0.2
    end

    @testset "XX point reachable via both LSM and XXZ delegation, agree" begin
        # Jx=Jy=1, Jz=0 hits the Jx=Jy branch (XXZ delegation) first; verify
        # the value matches XXZ1D(J=1, Δ=0) Energy(:per_site) directly.
        m_xx = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=0.0)
        f_dispatch = fetch(m_xx, GroundStateEnergyDensity(), Infinite())
        f_xxz = fetch(QAtlas.XXZ1D(; J=1.0, Δ=0.0), QAtlas.Energy{:per_site}(), Infinite())
        @test f_dispatch ≈ f_xxz atol=1e-12
    end

    @testset "Axial XXZ delegation: arbitrary Δ" begin
        for Δ in (-0.5, 0.0, 0.5, 0.9)
            m = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=Δ)
            f = fetch(m, GroundStateEnergyDensity(), Infinite())
            f_ref = fetch(
                QAtlas.XXZ1D(; J=1.0, Δ=Δ), QAtlas.Energy{:per_site}(), Infinite()
            )
            @test f ≈ f_ref atol=1e-12
        end
    end

    @testset "Generic XYZ (Jx≠Jy, Jz≠0) raises DomainError" begin
        m = HeisenbergXYZ(; Jx=2.0, Jy=1.0, Jz=0.5)
        @test_throws DomainError fetch(m, GroundStateEnergyDensity(), Infinite())
    end

    @testset "Axial XXZ with non-positive Jx raises DomainError" begin
        m = HeisenbergXYZ(; Jx=0.0, Jy=0.0, Jz=1.0)
        @test_throws DomainError fetch(m, GroundStateEnergyDensity(), Infinite())
        m_fm = HeisenbergXYZ(; Jx=-1.0, Jy=-1.0, Jz=0.5)
        @test_throws DomainError fetch(m_fm, GroundStateEnergyDensity(), Infinite())
    end

    @testset "Energy{:per_site} matches GroundStateEnergyDensity across regimes" begin
        using QAtlas: Energy
        for m in (
            HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=1.0),  # Heisenberg
            HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=0.5),  # XXZ
            HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=0.0),  # XX
            HeisenbergXYZ(; Jx=2.0, Jy=1.0, Jz=0.0),  # anisotropic XY
            HeisenbergXYZ(; Jx=3.0, Jy=0.5, Jz=0.0),  # extreme anisotropic XY
        )
            f_gsed = fetch(m, GroundStateEnergyDensity(), Infinite())
            f_ener = fetch(m, Energy{:per_site}(), Infinite())
            @test f_gsed == f_ener
        end
    end

    @testset "MassGap@Infinite -- critical axial + XY line" begin
        # Heisenberg, XXZ critical, XX point: all gapless
        for (Jx, Jy, Jz) in
            ((1.0, 1.0, 1.0), (1.0, 1.0, 0.5), (1.0, 1.0, -0.5), (1.0, 1.0, 0.0))
            m = HeisenbergXYZ(; Jx=Jx, Jy=Jy, Jz=Jz)
            @test fetch(m, QAtlas.MassGap(), Infinite()) == 0.0
        end
        # XY anisotropic line: gap = (1/4) |Jx - Jy|
        for (Jx, Jy) in ((2.0, 1.0), (3.0, 0.5), (1.0, 2.0))
            m = HeisenbergXYZ(; Jx=Jx, Jy=Jy, Jz=0.0)
            @test fetch(m, QAtlas.MassGap(), Infinite()) ≈ abs(Jx - Jy) / 4 atol=1e-12
        end
        # Massive AFM (Jz/Jx > 1): not yet implemented, DomainError
        m_massive = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=1.5)
        @test_throws DomainError fetch(m_massive, QAtlas.MassGap(), Infinite())
        # Generic XYZ: DomainError
        m_generic = HeisenbergXYZ(; Jx=2.0, Jy=1.0, Jz=0.5)
        @test_throws DomainError fetch(m_generic, QAtlas.MassGap(), Infinite())
    end

    @testset "CorrelationLength@Infinite -- critical + XY line" begin
        for (Jx, Jy, Jz) in ((1.0, 1.0, 1.0), (1.0, 1.0, 0.5), (1.0, 1.0, 0.0))
            m = HeisenbergXYZ(; Jx=Jx, Jy=Jy, Jz=Jz)
            @test fetch(m, QAtlas.CorrelationLength(), Infinite()) == Inf
        end
        for (Jx, Jy) in ((2.0, 1.0), (3.0, 0.5), (1.0, 2.0))
            m = HeisenbergXYZ(; Jx=Jx, Jy=Jy, Jz=0.0)
            xi_ref = 1 / asinh(abs(Jx - Jy) / (2 * sqrt(Jx * Jy)))
            @test fetch(m, QAtlas.CorrelationLength(), Infinite()) ≈ xi_ref atol=1e-12
        end
        m_strong = HeisenbergXYZ(; Jx=10.0, Jy=1.0, Jz=0.0)
        m_weak = HeisenbergXYZ(; Jx=2.0, Jy=1.0, Jz=0.0)
        @test fetch(m_strong, QAtlas.CorrelationLength(), Infinite()) <
            fetch(m_weak, QAtlas.CorrelationLength(), Infinite())
        m_massive = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=1.5)
        @test_throws DomainError fetch(m_massive, QAtlas.CorrelationLength(), Infinite())
        m_generic = HeisenbergXYZ(; Jx=2.0, Jy=1.0, Jz=0.5)
        @test_throws DomainError fetch(m_generic, QAtlas.CorrelationLength(), Infinite())
    end

    @testset "SpontaneousMagnetization@Infinite -- McCoy-Wu XY line + critical" begin
        for (Jx, Jy, Jz) in
            ((1.0, 1.0, 1.0), (1.0, 1.0, 0.5), (1.0, 1.0, -0.5), (1.0, 1.0, 0.0))
            m = HeisenbergXYZ(; Jx=Jx, Jy=Jy, Jz=Jz)
            @test fetch(m, QAtlas.SpontaneousMagnetization(), Infinite()) == 0.0
        end
        for (Jx, Jy) in ((2.0, 1.0), (3.0, 0.5), (10.0, 1.0), (1.0, 2.0))
            m = HeisenbergXYZ(; Jx=Jx, Jy=Jy, Jz=0.0)
            jmin, jmax = minmax(abs(Jx), abs(Jy))
            M_ref = (1 - (jmin / jmax)^2)^(1 / 8)
            @test fetch(m, QAtlas.SpontaneousMagnetization(), Infinite()) ≈ M_ref atol=1e-12
        end
        m_strong = HeisenbergXYZ(; Jx=100.0, Jy=1.0, Jz=0.0)
        @test fetch(m_strong, QAtlas.SpontaneousMagnetization(), Infinite()) > 0.999
        m_massive = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=1.5)
        @test_throws DomainError fetch(
            m_massive, QAtlas.SpontaneousMagnetization(), Infinite()
        )
        m_generic = HeisenbergXYZ(; Jx=2.0, Jy=1.0, Jz=0.5)
        @test_throws DomainError fetch(
            m_generic, QAtlas.SpontaneousMagnetization(), Infinite()
        )
    end

end
