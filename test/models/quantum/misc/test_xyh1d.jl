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
@testset "XYh1D — rejects Jx, Jy ≤ 0" begin
    @test_throws DomainError XYh1D(; Jx=0.0)
    @test_throws DomainError XYh1D(; Jx=-1.0)
    @test_throws DomainError XYh1D(; Jy=0.0)
    @test_throws DomainError XYh1D(; Jy=-1.0)
end
