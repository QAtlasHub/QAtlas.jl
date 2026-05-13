# ─────────────────────────────────────────────────────────────────────────────
# PpIp2DSC — 2-D p+ip chiral superconductor (Read-Green 2000), Phase 1.
#
# Verifies the two parameter-independent topological invariants of the
# weak-pairing phase:
#   * CentralCharge        c = 1/2   (chiral Majorana edge CFT)
#   * TopologicalInvariant C = 1     (first Chern number)
# and rejects parameters outside the weak-pairing branch.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "PpIp2DSC — edge CFT CentralCharge = 1/2 (Phase 1)" begin
    c = QAtlas.fetch(PpIp2DSC(), CentralCharge(), Infinite())
    @test c == 1 // 2
    @test c isa Rational
    # Parameter-independent within topological phase
    @test QAtlas.fetch(PpIp2DSC(; Δ₀=2.0, μ=3.0), CentralCharge(), Infinite()) == 1 // 2
    # Cross-check: chiral Majorana edge CFT = Ising boundary CFT = MinimalModel(4, 3)
    @test QAtlas.fetch(PpIp2DSC(), CentralCharge(), Infinite()) ==
        QAtlas.fetch(QAtlas.MinimalModel(4, 3), CentralCharge())
end

@testset "PpIp2DSC — TopologicalInvariant C = 1 (Phase 1)" begin
    @test QAtlas.fetch(PpIp2DSC(), TopologicalInvariant(), Infinite()) == 1
    @test QAtlas.fetch(PpIp2DSC(; Δ₀=5.0, μ=0.5), TopologicalInvariant(), Infinite()) == 1
end

@testset "PpIp2DSC — rejects Δ₀ ≤ 0 / μ ≤ 0 (Phase 1)" begin
    @test_throws DomainError PpIp2DSC(; Δ₀=0.0)
    @test_throws DomainError PpIp2DSC(; Δ₀=-1.0)
    @test_throws DomainError PpIp2DSC(; μ=0.0)
    @test_throws DomainError PpIp2DSC(; μ=-1.0)
end
