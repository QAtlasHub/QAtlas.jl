# ─────────────────────────────────────────────────────────────────────────────
# Tests for Compass1D — 1D alternating-bond compass chain (Phase 1).
#
# Coverage:
#   • Closed-form MassGap, Δ = 2 |J_x − J_y|, at several points.
#   • Symmetric point J_x = J_y: Δ = 0 (first-order QPT).
#   • DomainError on J_x ≤ 0 or J_y ≤ 0.
#
# References:
#   • Brzezicki–Dziarmaga–Oleś, PRB 75, 134415 (2007).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Compass1D — MassGap = 2|J_x − J_y| (Phase 1)" begin
    @test QAtlas.fetch(Compass1D(; J_x=1.0, J_y=0.5), MassGap(), Infinite()) == 1.0
    @test QAtlas.fetch(Compass1D(; J_x=2.0, J_y=0.5), MassGap(), Infinite()) == 3.0
    # Symmetric point: gap closes (first-order transition)
    @test QAtlas.fetch(Compass1D(; J_x=1.0, J_y=1.0), MassGap(), Infinite()) == 0.0
    @test QAtlas.fetch(Compass1D(; J_x=0.5, J_y=0.5), MassGap(), Infinite()) == 0.0
end

@testset "Compass1D — rejects J_x, J_y ≤ 0 (Phase 1)" begin
    @test_throws DomainError Compass1D(; J_x=0.0)
    @test_throws DomainError Compass1D(; J_x=-1.0)
    @test_throws DomainError Compass1D(; J_y=0.0)
end

@testset "Compass1D — swap symmetry: Δ(J_x, J_y) == Δ(J_y, J_x)" begin
    # True cross-check that the gap is |J_x − J_y| (symmetric in the
    # exchange of the two bond types), not an asymmetric subtraction.
    for (a, b) in ((1.0, 0.3), (0.7, 2.4), (1.5, 1.5), (0.01, 100.0))
        Δ_ab = QAtlas.fetch(Compass1D(; J_x=a, J_y=b), MassGap(), Infinite())
        Δ_ba = QAtlas.fetch(Compass1D(; J_x=b, J_y=a), MassGap(), Infinite())
        @test Δ_ab == Δ_ba
        @test Δ_ab == 2.0 * abs(a - b)
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "Compass1D — verification cards" begin
    for (jx, jy) in ((1.0, 0.5), (2.0, 0.5), (1.0, 1.0))
        verify(
            Compass1D(; J_x=jx, J_y=jy),
            MassGap(),
            Infinite();
            route=:second_closed_form,
            independent=2 * abs(jx - jy),
            agree_within=1e-10,
            refs=["1D compass model gap = 2|J_x - J_y|"],
        )
    end
end

# ── additional verification card (#381 batch) ─────────────────────────────
@testset "Compass1D — MassGap closed form (#381 batch)" begin
    # Brzezicki-Dziarmaga-Oles 2007: Δ = 2|J_x - J_y| (Jordan-Wigner dual of
    # the dimerised Kitaev chain). Gap closes at J_x = J_y (first-order QPT).
    # (1.0, 1.0) is included as a degenerate-gap regression guard: with
    # independent = 2*abs(0) = 0.0 exactly in IEEE 754, this also exercises
    # the hub at the critical point Δ = 0.
    for (Jx, Jy) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (0.7, 0.3), (2.0, 1.0))
        verify(
            Compass1D(; J_x=Jx, J_y=Jy),
            MassGap(),
            Infinite();
            route=:second_closed_form,
            independent=2 * abs(Jx - Jy),
            agree_within=1e-12,
            refs=["Brzezicki-Dziarmaga-Oles 2007 PRB 75 134415: Δ = 2|J_x - J_y| (JW-dual dimerised Kitaev chain)"],
        )
    end
end

