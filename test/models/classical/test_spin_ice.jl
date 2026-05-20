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

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "SpinIce — verification cards" begin
    # Pauling 1935 ice-rule residual entropy S = (1/2) log(3/2)
    verify(
        SpinIce(),
        ResidualEntropy(),
        Infinite();
        route=:second_closed_form,
        independent=0.5 * log(3 / 2),
        agree_within=1e-9,
        refs=["Pauling 1935: ice-rule residual entropy S = (1/2) log(3/2) ≈ 0.2027"],
    )
end
# ── additional verification cards (#381 batch 3) ─────────────────────────
@testset "SpinIce — Pauling ResidualEntropy (#381 batch 3)" begin
    # Pauling 1935: spin-ice T=0 residual entropy per spin is
    # (1/2) log(3/2) (half of the entropy per water molecule, since each
    # tetrahedron has 2 spins).
    verify(
        SpinIce(),
        ResidualEntropy(),
        Infinite();
        route=:second_closed_form,
        independent=log(3/2) / 2,
        agree_within=1e-12,
        refs=["Pauling 1935: spin-ice T=0 residual entropy = (1/2) log(3/2) per spin"],
    )
end

