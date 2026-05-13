# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: KitaevHeisenberg — K-only delegation to KitaevHoneycomb.
#
# Verifies:
#   * K-only isotropic point: MassGap = 0 (gapless Kitaev B-phase)
#   * J ≠ 0 raises DomainError (Phase 2 marker)
#   * Γ ≠ 0 raises DomainError (Phase 2 marker)
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "KitaevHeisenberg — K-only isotropic delegates to gapless Kitaev B" begin
    m = KitaevHeisenberg(; K=1.0, J=0.0, Γ=0.0)
    Δ = QAtlas.fetch(m, MassGap(), Infinite())
    @test Δ == 0.0
end

@testset "KitaevHeisenberg — K-only K-scaling" begin
    for K in (0.5, 1.0, 2.5)
        m = KitaevHeisenberg(; K=K, J=0.0, Γ=0.0)
        # Isotropic K-only honeycomb is in the gapless B phase: Δ = 0.
        @test QAtlas.fetch(m, MassGap(), Infinite()) == 0.0
    end
end

@testset "KitaevHeisenberg — DomainError on J ≠ 0" begin
    @test_throws DomainError QAtlas.fetch(
        KitaevHeisenberg(; K=1.0, J=0.1, Γ=0.0), MassGap(), Infinite()
    )
end

@testset "KitaevHeisenberg — DomainError on Γ ≠ 0" begin
    @test_throws DomainError QAtlas.fetch(
        KitaevHeisenberg(; K=1.0, J=0.0, Γ=0.1), MassGap(), Infinite()
    )
end

@testset "KitaevHeisenberg — DomainError on both J and Γ" begin
    @test_throws DomainError QAtlas.fetch(
        KitaevHeisenberg(; K=1.0, J=0.1, Γ=0.1), MassGap(), Infinite()
    )
end
