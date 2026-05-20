# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: TodaLattice — classical Toda gapless phonon.
#
# Verifies:
#   * MassGap = 0 for all (a, b) > 0  (acoustic phonon at k = 0)
#   * (a, b) defaults yield same result as explicit construction
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TodaLattice — gapless acoustic phonon MassGap = 0" begin
    @test QAtlas.fetch(TodaLattice(), MassGap(), Infinite()) == 0.0
    @test QAtlas.fetch(TodaLattice(; a=1.0, b=1.0), MassGap(), Infinite()) == 0.0
    for a in (0.5, 2.0, 3.7), b in (0.5, 1.0, 2.5)
        @test QAtlas.fetch(TodaLattice(; a=a, b=b), MassGap(), Infinite()) == 0.0
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "TodaLattice — verification cards" begin
    # The classical Toda chain has a gapless acoustic phonon branch:
    # MassGap = 0 for all (a, b > 0) (integrable dispersion at q -> 0).
    for (a, b) in ((1.0, 1.0), (2.0, 0.5), (0.7, 1.3))
        verify(
            TodaLattice(; a=a, b=b),
            MassGap(),
            Infinite();
            route=:second_closed_form,
            independent=0.0,
            agree_within=1e-12,
            refs=["Toda chain acoustic branch is gapless: MassGap = 0"],
        )
    end
end
# ── additional verification cards (#381 batch 3) ─────────────────────────
@testset "TodaLattice — MassGap gapless (#381 batch 3)" begin
    # Classical Toda lattice is integrable (Flaschka 1974; Henon 1974) with
    # acoustic phonons at low momentum ⇒ no gap, Δ = 0.
    verify(
        TodaLattice(),
        MassGap(),
        Infinite();
        route=:second_closed_form,
        independent=0.0,
        agree_within=1e-12,
        refs=["Flaschka 1974 / Henon 1974: integrable Toda lattice with acoustic phonons ⇒ MassGap = 0"],
    )
end

