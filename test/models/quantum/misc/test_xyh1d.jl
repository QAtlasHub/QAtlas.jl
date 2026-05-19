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

# ─────────────────────────────────────────────────────────────────────────────
# Energy{:per_site} — closed-form XX limit (Jx = Jy), any h.  LSM 1961.
# ─────────────────────────────────────────────────────────────────────────────

@testset "XYh1D — XX limit Energy{:per_site} (Phase 1)" begin
    # h = 0, J = 1: E/N = -4/π  (Lieb-Schultz-Mattis 1961)
    e0 = QAtlas.fetch(XYh1D(), Energy{:per_site}(), Infinite())
    @test isapprox(e0, -4 / π; atol=1e-12)
    @test isapprox(e0, -1.2732395447351628; atol=1e-12)

    # J linearity at h = 0
    for J in (0.5, 1.0, 2.0, 3.5)
        e = QAtlas.fetch(XYh1D(; Jx=J, Jy=J, h=0.0), Energy{:per_site}(), Infinite())
        @test isapprox(e, -4J / π; atol=1e-12)
    end

    # Fully polarised (|h| ≥ 2J): E/N = -|h|
    @test QAtlas.fetch(XYh1D(; h=2.0), Energy{:per_site}(), Infinite()) ≈ -2.0
    @test QAtlas.fetch(XYh1D(; h=3.0), Energy{:per_site}(), Infinite()) ≈ -3.0
    @test QAtlas.fetch(XYh1D(; h=-2.0), Energy{:per_site}(), Infinite()) ≈ -2.0
    @test QAtlas.fetch(XYh1D(; h=-5.0), Energy{:per_site}(), Infinite()) ≈ -5.0

    # Continuity at h = 2J: formula matches polarised value
    e_just_below = QAtlas.fetch(XYh1D(; h=2.0 - 1e-9), Energy{:per_site}(), Infinite())
    @test isapprox(e_just_below, -2.0; atol=1e-4)

    # Inside band: closed form against direct evaluation
    let J = 1.0, h = 1.0
        x = h / (2J)
        ref = -h + (2h / π) * acos(x) - (4J / π) * sqrt(1 - x^2)
        @test QAtlas.fetch(XYh1D(; h=h), Energy{:per_site}(), Infinite()) ≈ ref
    end
end

@testset "XYh1D — Energy{:per_site} anisotropic throws (Phase 2)" begin
    @test_throws DomainError QAtlas.fetch(
        XYh1D(; Jx=1.0, Jy=0.5, h=0.0), Energy{:per_site}(), Infinite()
    )
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "XYh1D — verification cards" begin
    verify(
        XYh1D(),
        Energy(:per_site),
        Infinite();
        route=:second_closed_form,
        independent=-4 / pi,
        agree_within=1e-9,
        refs=["Lieb-Schultz-Mattis 1961: XX chain e0 = -4/pi (Pauli σ convention)"],
    )
    verify(
        XYh1D(),
        MassGap(),
        Infinite();
        route=:second_closed_form,
        independent=0.0,
        agree_within=1e-10,
        refs=["XX limit h=0: gapless"],
    )
    verify(
        XYh1D(; h=3.0),
        MassGap(),
        Infinite();
        route=:second_closed_form,
        independent=2.0,
        agree_within=1e-9,
        refs=["Polarized |h|>2J: gap = 2(|h| - 2J)"],
    )
end
