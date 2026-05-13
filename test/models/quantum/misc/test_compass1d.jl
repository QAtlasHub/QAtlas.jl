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
