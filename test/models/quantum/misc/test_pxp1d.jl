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
    # Linear in Ω
    @test isapprox(
        QAtlas.fetch(PXP1D(; Ω=2.5), Energy{:per_site}(), Infinite()),
        2.5 * (-0.6516);
        atol=2.5e-4,
    )
end

@testset "PXP1D — rejects Ω ≤ 0 (Phase 1)" begin
    @test_throws DomainError PXP1D(; Ω=0.0)
    @test_throws DomainError PXP1D(; Ω=-1.0)
    m = PXP1D(; Ω=1.0)
    @test_throws DomainError QAtlas.fetch(m, Energy{:per_site}(), Infinite(); Ω=0.0)
end
