# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: TricriticalPotts3 — c = 6/7 via MinimalModel(6,7).
#
# Verifies:
#   * Central charge exactly 6//7 (Rational, machine-precision agreement).
#   * Result equals MinimalModel(7, 6)'s c (delegation invariant).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TricriticalPotts3 — c = 6/7 exact" begin
    c = QAtlas.fetch(TricriticalPotts3(), CentralCharge(), Infinite())
    @test c == 6 // 7
end

@testset "TricriticalPotts3 — equals MinimalModel(7, 6)" begin
    c_tp = QAtlas.fetch(TricriticalPotts3(), CentralCharge(), Infinite())
    c_mm = QAtlas.fetch(QAtlas.MinimalModel(7, 6), CentralCharge(), Infinite())
    @test c_tp == c_mm
end
