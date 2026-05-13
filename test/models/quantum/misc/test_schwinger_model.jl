# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: SchwingerModel — massless Schwinger m_γ = e/√π.
#
# Verifies:
#   * m_γ = e/√π exactly at m = 0 for several e values
#   * Scaling: m_γ(2e) = 2 m_γ(e)
#   * DomainError on e ≤ 0
#   * DomainError on m ≠ 0 (massive case is Phase 2)
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "SchwingerModel — massless m_γ = e/√π" begin
    for e in (0.5, 1.0, 2.0, 3.7)
        mγ = QAtlas.fetch(SchwingerModel(; e=e, m=0.0), MassGap(), Infinite())
        @test mγ ≈ e / sqrt(π) atol = 1e-14
    end
end

@testset "SchwingerModel — linear e-scaling" begin
    m1 = QAtlas.fetch(SchwingerModel(; e=1.0), MassGap(), Infinite())
    m2 = QAtlas.fetch(SchwingerModel(; e=2.0), MassGap(), Infinite())
    @test m2 ≈ 2 * m1 atol = 1e-14
end

@testset "SchwingerModel — DomainError on e ≤ 0" begin
    @test_throws DomainError QAtlas.fetch(SchwingerModel(; e=0.0), MassGap(), Infinite())
    @test_throws DomainError QAtlas.fetch(SchwingerModel(; e=-1.0), MassGap(), Infinite())
end

@testset "SchwingerModel — DomainError on m ≠ 0 (Phase 2)" begin
    @test_throws DomainError QAtlas.fetch(
        SchwingerModel(; e=1.0, m=0.1), MassGap(), Infinite()
    )
end
