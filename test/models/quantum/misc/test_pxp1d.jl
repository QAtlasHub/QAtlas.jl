# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: PXP1D — Rydberg-blockade chain (Phase 1).
#
# Verifies:
#   * Energy{:per_site} at Infinite returns the hardcoded Surace 2020
#     MPS thermodynamic-limit reference  e_0 / Ω ≈ -0.6516.
#   * Linearity in the coupling Ω.
#   * DomainError on non-positive Ω (constructor + fetch override).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "PXP1D — Energy{:per_site} reference value (Phase 1)" begin
    e0 = QAtlas.fetch(PXP1D(), Energy{:per_site}(), Infinite())
    @test e0 ≈ -0.6516
    @test e0 < 0
    # Linear in Ω
    @test QAtlas.fetch(PXP1D(; Ω=2.5), Energy{:per_site}(), Infinite()) ≈ 2.5 * (-0.6516)
end

@testset "PXP1D — rejects Ω ≤ 0 (Phase 1)" begin
    @test_throws DomainError PXP1D(; Ω=0.0)
    @test_throws DomainError PXP1D(; Ω=-1.0)
    m = PXP1D(; Ω=1.0)
    @test_throws DomainError QAtlas.fetch(m, Energy{:per_site}(), Infinite(); Ω=0.0)
end
