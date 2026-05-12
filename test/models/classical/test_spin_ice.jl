# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: SpinIce — Pauling 1935 residual entropy.
#
# Verifies:
#   * S/N matches the Pauling closed form (1/2) log(3/2) to 1e-14
#   * Pauling value is below the Nagle 1966 numerical S/N ≈ 0.20479
#     by 0.5 %–1 % — pinned as a sanity-check inequality
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "SpinIce — Pauling residual entropy" begin
    S_per_spin = QAtlas.fetch(SpinIce(), ResidualEntropy(), Infinite())
    @test S_per_spin ≈ 0.5 * log(3 / 2) atol = 1e-14
    @test S_per_spin ≈ 0.2027325540540822 atol = 1e-12
    # Nagle 1966 numerical result for the pyrochlore ice rule:
    # S/N ≈ 0.20479; Pauling underestimates by 0.5 %.
    @test S_per_spin < 0.20479
    @test S_per_spin > 0.20         # comfortably above the Pauling −1 % band
end
