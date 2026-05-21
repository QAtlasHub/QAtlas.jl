# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: PXP1D — Rydberg-blockade chain (Phase 1).
#
# Verifies:
#   * Energy{:per_site} at Infinite returns the hardcoded DMRG/ED
#     thermodynamic-limit reference  e_0 / Ω ≈ -0.6516 (PXP-scar
#     literature; see Turner 2018, Lin-Motrunich 2019, Iadecola 2019).
#   * Linearity in the coupling Ω.
#   * DomainError on non-positive Ω (constructor + fetch override).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "PXP1D — Energy{:per_site} reference value (Phase 1)" begin
    e0 = QAtlas.fetch(PXP1D(), Energy{:per_site}(), Infinite())
    @test isapprox(e0, -0.6516; atol=1e-4)  # literature precision -0.6516(2)
    @test e0 < 0
    @test e0 isa Float64
    # Linear in Ω
    @test isapprox(
        QAtlas.fetch(PXP1D(; Ω=2.5), Energy{:per_site}(), Infinite()),
        2.5 * (-0.6516);
        atol=2.5e-4,
    )
    # Sign preserved across scales (Ω > 0 ⇒ e_0 < 0)
    @test QAtlas.fetch(PXP1D(; Ω=10.0), Energy{:per_site}(), Infinite()) < 0
    # Exact linearity ratio (no hidden additive constant)
    e1 = QAtlas.fetch(PXP1D(; Ω=1.0), Energy{:per_site}(), Infinite())
    e7 = QAtlas.fetch(PXP1D(; Ω=7.0), Energy{:per_site}(), Infinite())
    @test isapprox(e7 / e1, 7.0; atol=1e-12)
end

@testset "PXP1D — rejects Ω ≤ 0 (Phase 1)" begin
    @test_throws DomainError PXP1D(; Ω=0.0)
    @test_throws DomainError PXP1D(; Ω=-1.0)
    m = PXP1D(; Ω=1.0)
    @test_throws DomainError QAtlas.fetch(m, Energy{:per_site}(), Infinite(); Ω=0.0)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "PXP1D — verification cards" begin
    verify(
        PXP1D(),
        Energy(:per_site),
        Infinite();
        route=:literature_value,
        independent=-0.6516,
        agree_within=5e-3,
        refs=["PXP scar literature (Turner 2018 / Lin-Motrunich 2019): e0 ~ -0.6516"],
    )
end

# ── additional verification cards (#381 batch 4) ─────────────────────────
@testset "PXP1D — Energy DMRG reference (#381 batch 4)" begin
    # PXP model GS energy density e_0 ≈ -0.6516 (Lin-Motrunich 2019 DMRG,
    # PRL 122, 173401; matches src tabulated value).
    verify(
        PXP1D(),
        Energy(:per_site),
        Infinite();
        route=:literature_value,
        independent=-0.6516,
        agree_within=1e-3,
        refs=["Lin-Motrunich 2019 PRL 122 173401: PXP DMRG GS energy density ≈ -0.6516"],
    )
end
