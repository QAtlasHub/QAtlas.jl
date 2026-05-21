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
    for Jp in (1.0, 2.0, 0.5), α in (0.0, 0.1, 0.3, 0.6, 0.675)
        J = α * Jp
        m = QAtlas.ShastrySutherland(; J=J, Jp=Jp)
        E0 = QAtlas.fetch(m, Energy(:per_site), Infinite())
        @test E0 ≈ -3 * Jp / 8 atol = 1e-14
    end
end

@testset "ShastrySutherland — independence of J inside dimer window" begin
    Jp = 1.0
    vals = [
        QAtlas.fetch(QAtlas.ShastrySutherland(; J=J, Jp=Jp), Energy(:per_site), Infinite())
        for J in (0.0, 0.2, 0.4, 0.65)
    ]
    @test all(v -> isapprox(v, -3 / 8, atol=1e-14), vals)
end

@testset "ShastrySutherland — DomainError outside dimer phase" begin
    m_large_alpha = QAtlas.ShastrySutherland(; J=0.7, Jp=1.0)   # α = 0.7 > α_c
    @test_throws DomainError QAtlas.fetch(m_large_alpha, Energy(:per_site), Infinite())
end

@testset "ShastrySutherland — DomainError on non-AF Jp" begin
    @test_throws DomainError QAtlas.fetch(
        QAtlas.ShastrySutherland(; J=0.0, Jp=0.0), Energy(:per_site), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        QAtlas.ShastrySutherland(; J=0.0, Jp=-1.0), Energy(:per_site), Infinite()
    )
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "ShastrySutherland — verification cards" begin
    for (J, Jp) in ((0.1, 1.0), (0.6, 2.0), (0.3, 1.5))
        verify(
            QAtlas.ShastrySutherland(; J=J, Jp=Jp),
            Energy(:per_site),
            Infinite();
            route=:second_closed_form,
            independent=-3 * Jp / 8,
            agree_within=1e-10,
            refs=["Shastry-Sutherland 1981: exact dimer e0 = -3 Jp / 8 (dimer phase)"],
        )
    end
end

# ── additional verification card (#381 batch) ─────────────────────────────
@testset "ShastrySutherland — Energy/Infinite dimer-phase card (#381 batch)" begin
    # Shastry-Sutherland 1981 exact dimer phase: for α = J/Jp ≤ α_c ≈ 0.675
    # (Koga-Kawakami 2000) the GS is the exact product of singlets on the
    # diagonal J' bonds; orthogonality cancels the J nearest-neighbour bond
    # contribution, giving e₀ = -3 J' / 8 independent of J.
    # agree_within=1e-14 is tighter than the file’s pre-existing 1e-10 floor:
    # both card (-3*Jp/8) and hub evaluate the same single multiply+divide, so the
    # closed-form path is bit-identical for all sweep points in IEEE 754.
    # Sweep includes (0.674, 1.0) just below the dimer-phase boundary α_c ≈ 0.675
    # (Koga-Kawakami 2000) as a regression guard for the phase-boundary logic.
    for (J, Jp) in
        ((0.0, 1.0), (0.3, 1.0), (0.674, 1.0), (0.0, 2.0), (0.5, 2.0), (1.0, 2.0))
        verify(
            QAtlas.ShastrySutherland(; J=J, Jp=Jp),
            Energy(:per_site),
            Infinite();
            route=:second_closed_form,
            independent=-3 * Jp / 8,
            agree_within=1e-14,
            refs=[
                "Shastry-Sutherland 1981 / Koga-Kawakami 2000: exact dimer phase e₀ = -3J'/8 for α ≤ α_c ≈ 0.675",
            ],
        )
    end
end
