# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: XYh1D — isotropic XX limit MassGap = 2·max(0, |h| − 2J)
# (Lieb-Schultz-Mattis 1961; Pfeuty 1970).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "XYh1D — isotropic XX limit MassGap (Phase 1)" begin
    # h = 0: gapless XX chain
    @test QAtlas.fetch(XYh1D(), MassGap(), Infinite()) == 0.0
    # Inside band (|h| < 2J): gapless
    @test QAtlas.fetch(XYh1D(; h=1.0), MassGap(), Infinite()) == 0.0
    # At the critical h = 2J: gap closes (Lifshitz / BKT-like point)
    @test QAtlas.fetch(XYh1D(; h=2.0), MassGap(), Infinite()) == 0.0
    # Polarised (|h| > 2J): finite gap = 2(|h| - 2J)
    @test QAtlas.fetch(XYh1D(; h=3.0), MassGap(), Infinite()) == 2.0
    @test QAtlas.fetch(XYh1D(; h=-3.0), MassGap(), Infinite()) == 2.0  # depends on |h|
    @test QAtlas.fetch(XYh1D(; h=5.0), MassGap(), Infinite()) == 6.0
    # Different J (still isotropic Jx = Jy)
    @test QAtlas.fetch(XYh1D(; Jx=0.5, Jy=0.5, h=2.0), MassGap(), Infinite()) == 2.0
end

@testset "XYh1D — anisotropic case throws DomainError (deferred to Phase 2)" begin
    @test_throws DomainError QAtlas.fetch(
        XYh1D(; Jx=1.0, Jy=0.5, h=0.5), MassGap(), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        XYh1D(; Jx=2.0, Jy=1.0, h=2.0), MassGap(), Infinite()
    )
end

@testset "XYh1D — rejects Jx, Jy ≤ 0 (Phase 1)" begin
    @test_throws DomainError XYh1D(; Jx=0.0)
    @test_throws DomainError XYh1D(; Jx=-1.0)
    @test_throws DomainError XYh1D(; Jy=0.0)
    @test_throws DomainError XYh1D(; Jy=-1.0)
end
