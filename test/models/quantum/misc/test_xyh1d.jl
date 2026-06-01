# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: XYh1D — general anisotropic XY model exact solutions
# (Lieb-Schultz-Mattis 1961; Pfeuty 1970).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

if !isdefined(Main, :verify)
    include("../../../util/verify.jl")
end

@testset "XYh1D — MassGap (Infinite and OBC)" begin
    # h = 0: gapless XX chain
    @test QAtlas.fetch(XYh1D(), MassGap(), Infinite()) == 0.0
    # Inside band (|h| < 2J): gapless
    @test QAtlas.fetch(XYh1D(; h=1.0), MassGap(), Infinite()) == 0.0
    # At the critical h = 2J: gap closes
    @test QAtlas.fetch(XYh1D(; h=2.0), MassGap(), Infinite()) == 0.0
    # Polarised (|h| > 2J): finite gap = 2(|h| - 2J)
    @test QAtlas.fetch(XYh1D(; h=3.0), MassGap(), Infinite()) == 2.0
    @test QAtlas.fetch(XYh1D(; h=-3.0), MassGap(), Infinite()) == 2.0
    @test QAtlas.fetch(XYh1D(; h=5.0), MassGap(), Infinite()) == 6.0

    # Anisotropic Jx ≠ Jy (MassGap should be finite for h < Jx+Jy if Jx ≠ Jy)
    # e.g., Jx=1.0, Jy=0.5, h=0.0: MassGap = 2 * |Jx - Jy| = 2.0 * 0.5 = 1.0
    @test QAtlas.fetch(XYh1D(; Jx=1.0, Jy=0.5, h=0.0), MassGap(), Infinite()) ≈ 1.0

    # OBC MassGap matches the smallest positive eigenvalue of the BdG spectrum
    let m = XYh1D(; Jx=1.0, Jy=0.5, h=0.5)
        gap_inf = QAtlas.fetch(m, MassGap(), Infinite())
        gap_obc = QAtlas.fetch(m, MassGap(), OBC(200))
        @test isapprox(gap_obc, gap_inf; atol=1e-3)
    end
end

@testset "XYh1D — rejects Jx, Jy ≤ 0" begin
    @test_throws DomainError XYh1D(; Jx=0.0)
    @test_throws DomainError XYh1D(; Jx=-1.0)
    @test_throws DomainError XYh1D(; Jy=0.0)
    @test_throws DomainError XYh1D(; Jy=-1.0)
end

@testset "XYh1D — Energy{:per_site} (Infinite and OBC)" begin
    # h = 0, J = 1: E/N = -4/π  (Lieb-Schultz-Mattis 1961)
    e0 = QAtlas.fetch(XYh1D(), Energy{:per_site}(), Infinite())
    @test isapprox(e0, -4 / π; atol=1e-12)

    # J linearity at h = 0
    for J in (0.5, 1.0, 2.0, 3.5)
        e = QAtlas.fetch(XYh1D(; Jx=J, Jy=J, h=0.0), Energy{:per_site}(), Infinite())
        @test isapprox(e, -4J / π; atol=1e-12)
    end
end

@testset "XYh1D — Finite Temperature potentials" begin
    # Check that OBC free energy approaches the infinite thermodynamic limit
    let m = XYh1D(; Jx=1.2, Jy=0.8, h=0.6), β = 2.0
        f_inf = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
        f_obc = QAtlas.fetch(m, FreeEnergy(), OBC(200); beta=β)
        @test isapprox(f_obc, f_inf; atol=1e-2)

        s_inf = QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=β)
        s_obc = QAtlas.fetch(m, ThermalEntropy(), OBC(200); beta=β)
        @test isapprox(s_obc, s_inf; atol=1e-2)

        c_inf = QAtlas.fetch(m, SpecificHeat(), Infinite(); beta=β)
        c_obc = QAtlas.fetch(m, SpecificHeat(), OBC(200); beta=β)
        @test isapprox(c_obc, c_inf; atol=1e-2)
    end
end
