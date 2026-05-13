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
