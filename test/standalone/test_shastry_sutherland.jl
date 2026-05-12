# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: ShastrySutherland — exact dimer ground state.
#
# Verifies:
#   * E0/N = -3 J'/8 exactly across the dimer phase J/J' ≤ α_c ≈ 0.675
#   * Result is independent of J in the dimer window (Shastry-Sutherland
#     1981 / Koga-Kawakami 2000)
#   * J'-scaling: doubling J' doubles -E0
#   * DomainError outside the dimer phase
#   * DomainError on Jp ≤ 0 (non-AF dimer bond)
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "ShastrySutherland — exact dimer E0/N = -3J'/8" begin
    for Jp in (1.0, 2.0, 0.5), J in (0.0, 0.1, 0.3, 0.6, 0.675)
        m = ShastrySutherland(; J=J, Jp=Jp)
        E0 = QAtlas.fetch(m, Energy(:per_site), Infinite())
        @test E0 ≈ -3 * Jp / 8 atol = 1e-14
    end
end

@testset "ShastrySutherland — independence of J inside dimer window" begin
    Jp = 1.0
    vals = [
        QAtlas.fetch(ShastrySutherland(; J=J, Jp=Jp), Energy(:per_site), Infinite()) for
        J in (0.0, 0.2, 0.4, 0.65)
    ]
    @test all(v -> isapprox(v, -3 / 8, atol=1e-14), vals)
end

@testset "ShastrySutherland — DomainError outside dimer phase" begin
    m_large_alpha = ShastrySutherland(; J=0.7, Jp=1.0)   # α = 0.7 > α_c
    @test_throws DomainError QAtlas.fetch(m_large_alpha, Energy(:per_site), Infinite())
end

@testset "ShastrySutherland — DomainError on non-AF Jp" begin
    @test_throws DomainError QAtlas.fetch(
        ShastrySutherland(; J=0.0, Jp=0.0), Energy(:per_site), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        ShastrySutherland(; J=0.0, Jp=-1.0), Energy(:per_site), Infinite()
    )
end
