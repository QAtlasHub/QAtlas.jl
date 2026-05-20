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

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "PpIp2DSC — verification cards" begin
    verify(
        PpIp2DSC(),
        CentralCharge(),
        Infinite();
        route=:second_closed_form,
        independent=1 // 2,
        agree_within=1e-10,
        refs=["Chiral p+ip edge: single chiral Majorana CFT c = 1/2"],
    )
    verify(
        PpIp2DSC(),
        TopologicalInvariant(),
        Infinite();
        route=:second_closed_form,
        independent=1,
        agree_within=1e-10,
        refs=["p+ip weak-pairing phase: first Chern number = 1"],
    )
end
# ── additional verification cards (#381 batch 3) ─────────────────────────
@testset "PpIp2DSC — CentralCharge (#381 batch 3)" begin
    # Chiral p+ip 2D superconductor has Majorana edge mode ⇒ boundary CFT
    # is the chiral Ising/free-Majorana with c = 1/2 (Read-Green 2000;
    # Kitaev 2003).
    verify(
        PpIp2DSC(),
        CentralCharge(),
        Infinite();
        route=:second_closed_form,
        independent=1//2,
        agree_within=0,
        refs=["Read-Green 2000; Kitaev 2003: chiral p+ip SC has Majorana edge ⇒ c = 1/2"],
    )
end
