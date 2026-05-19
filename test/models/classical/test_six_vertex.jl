# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: SixVertex — Lieb 1967 ice rule + KDP/F-model phases.
#
# Verifies (Phase 1 + Phase 2 of issue #163):
#
#   * Square ice (a = b = c = 1) ResidualEntropy = (3/2) log(4/3)
#     ≈ 0.4315231087… at atol 1e-14 (Lieb 1967a, exact closed form).
#
#   * Phase classification via Δ = (a² + b² − c²) / (2 a b):
#       Δ > 1   → :ferroelectric
#       |Δ| ≤ 1 → :disordered
#       Δ < −1  → :antiferroelectric
#
#   * Ferroelectric phase (Δ > 1):
#       ResidualEntropy = 0
#       FreeEnergy      = −log max(a, b)
#     Verified at the KDP point (a = 2, b = c = 1) and at a symmetric
#     b-dominated point (a = 1, b = 3, c = 1).
#
#   * Antiferroelectric phase (Δ < −1, e.g. f-model c = 3):
#       FreeEnergy currently throws ArgumentError (deferred to phase 3
#       of issue #163, Lieb 1967b elliptic-function form).
#     This is asserted explicitly so the deferral is visible.
#
#   * Generic disordered point off the square-ice diagonal
#     (a = b = 1, c = 0.5, Δ = 0.875): FreeEnergy currently throws
#     ArgumentError (deferred to phase 2 of issue #163).  Asserted
#     explicitly so the scope boundary is documented.
#
#   * Square-ice FreeEnergy = −(3/2) log(4/3) closed form
#     (Lieb 1967a, dual of the residual entropy at zero temperature).
#
# No Lattice2D dependency.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "SixVertex — square ice residual entropy (Lieb 1967a)" begin
    m = QAtlas.square_ice()
    @test m isa SixVertex
    @test m.a == 1.0 && m.b == 1.0 && m.c == 1.0

    S = QAtlas.fetch(m, ResidualEntropy(), Infinite())
    S_lieb = (3 / 2) * log(4 / 3)
    @test S ≈ S_lieb atol = 1e-14
    # The exact value of (3/2) log(4/3) is 0.4315231086776713 (the
    # ≈ 0.4314 quoted in some surveys is a rounding of the truncated
    # 4-digit form; cross-check at the wider tolerance the issue
    # description used).
    @test S ≈ 0.4314 atol = 1e-3
    @test S ≈ 0.43152310867767 atol = 1e-12
end

@testset "SixVertex — phase classification" begin
    # Square ice: Δ = 1/2 ⇒ disordered.
    @test QAtlas._six_vertex_delta(1.0, 1.0, 1.0) ≈ 0.5
    @test QAtlas._six_vertex_phase(1.0, 1.0, 1.0) === :disordered

    # KDP critical point a = 2: Δ = a/2 = 1 (boundary, classified
    # disordered).  Note: with b = c = 1, Δ = a²/(2a) = a/2, so the
    # KDP transition is at a = 2, *not* √2.
    @test QAtlas._six_vertex_delta(2.0, 1.0, 1.0) ≈ 1.0 atol = 1e-14
    @test QAtlas._six_vertex_phase(2.0, 1.0, 1.0) === :disordered

    # KDP frozen FE: a = 3, b = c = 1 ⇒ Δ = 3/2 > 1.
    @test QAtlas._six_vertex_delta(3.0, 1.0, 1.0) > 1.0
    @test QAtlas._six_vertex_phase(3.0, 1.0, 1.0) === :ferroelectric

    # f-model critical point c = 2: Δ = (1 + 1 − 4)/2 = −1.
    @test QAtlas._six_vertex_delta(1.0, 1.0, 2.0) ≈ -1.0
    @test QAtlas._six_vertex_phase(1.0, 1.0, 2.0) === :disordered

    # f-model AFE: c = 3 ⇒ Δ = (1 + 1 − 9)/2 = −7/2.
    @test QAtlas._six_vertex_delta(1.0, 1.0, 3.0) ≈ -7 / 2
    @test QAtlas._six_vertex_phase(1.0, 1.0, 3.0) === :antiferroelectric
end

@testset "SixVertex — ferroelectric phase (Lieb 1967c)" begin
    # KDP frozen FE: a > 2, b = c = 1.  S_residual = 0, f = −log a.
    m_kdp = QAtlas.kdp_model(3.0)
    @test QAtlas._six_vertex_phase(m_kdp.a, m_kdp.b, m_kdp.c) === :ferroelectric

    @test QAtlas.fetch(m_kdp, ResidualEntropy(), Infinite()) == 0.0
    @test QAtlas.fetch(m_kdp, FreeEnergy(), Infinite()) ≈ -log(3.0) atol = 1e-14

    # Symmetric b-dominated FE: a = 1, b = 3, c = 1.  Δ = (1+9-1)/(2·3) = 9/6 = 1.5 > 1.
    m_b = SixVertex(; a=1.0, b=3.0, c=1.0)
    @test QAtlas._six_vertex_phase(m_b.a, m_b.b, m_b.c) === :ferroelectric
    @test QAtlas.fetch(m_b, ResidualEntropy(), Infinite()) == 0.0
    @test QAtlas.fetch(m_b, FreeEnergy(), Infinite()) ≈ -log(3.0) atol = 1e-14

    # Sanity: kdp_model(a) only enters FE for a > 2.
    m_below = QAtlas.kdp_model(1.5)        # a = 1.5 < 2, still disordered.
    @test QAtlas._six_vertex_phase(m_below.a, m_below.b, m_below.c) === :disordered
end

@testset "SixVertex — antiferroelectric phase (deferred — Lieb 1967b)" begin
    # f-model with c = 3 ⇒ Δ = −7/2 < −1, AFE phase.
    m_afe = QAtlas.f_model(3.0)
    @test QAtlas._six_vertex_phase(m_afe.a, m_afe.b, m_afe.c) === :antiferroelectric

    # Phase 3 of issue #163: AFE elliptic-function free energy is not
    # implemented in this commit.  The fetch must throw an informative
    # ArgumentError so deferred coverage is visible.
    @test_throws ArgumentError QAtlas.fetch(m_afe, FreeEnergy(), Infinite())

    # Likewise, generic-disordered + AFE residual entropy is deferred.
    @test_throws ArgumentError QAtlas.fetch(m_afe, ResidualEntropy(), Infinite())
end

@testset "SixVertex — disordered free energy (square-ice closed form)" begin
    # Square-ice point a = b = c = 1: Lieb 1967a gives
    # f = -(3/2) log(4/3) (the negative of the residual entropy).
    m_si = QAtlas.square_ice()
    f_si = QAtlas.fetch(m_si, FreeEnergy(), Infinite())
    @test f_si ≈ -(3 / 2) * log(4 / 3) atol = 1e-14
    # At zero temperature f and -S coincide for the ice-rule manifold.
    @test f_si ≈ -QAtlas.fetch(m_si, ResidualEntropy(), Infinite()) atol = 1e-14
end

@testset "SixVertex — disordered free energy off-diagonal (deferred — phase 2)" begin
    # Generic disordered point (a = b = 1, c = 0.5):
    #   Δ = (1 + 1 − 0.25) / 2 = 0.875 ∈ (-1, 1).
    # The Lieb / Sutherland 1967 trigonometric integral covering the
    # full disordered phase is deferred to a follow-up commit (issue
    # #163 phase 2); the fetch must throw an informative ArgumentError.
    m_d = SixVertex(; a=1.0, b=1.0, c=0.5)
    @test QAtlas._six_vertex_delta(1.0, 1.0, 0.5) ≈ 0.875
    @test QAtlas._six_vertex_phase(m_d.a, m_d.b, m_d.c) === :disordered
    @test_throws ArgumentError QAtlas.fetch(m_d, FreeEnergy(), Infinite())

    # f-model boundary c = 2 (Δ = −1, disordered branch but a = b ≠ c):
    # also deferred.
    m_crit = QAtlas.f_model(2.0)
    @test QAtlas._six_vertex_phase(m_crit.a, m_crit.b, m_crit.c) === :disordered
    @test_throws ArgumentError QAtlas.fetch(m_crit, FreeEnergy(), Infinite())
end

@testset "SixVertex — constructor argument validation" begin
    @test_throws ArgumentError SixVertex(; a=0.0)
    @test_throws ArgumentError SixVertex(; b=-1.0)
    @test_throws ArgumentError SixVertex(; c=-0.5)
    # Default constructor lands on square ice.
    @test SixVertex() === SixVertex(1.0, 1.0, 1.0)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "SixVertex — verification cards" begin
    # Square ice (a=b=c=1): Lieb 1967 residual entropy S = (3/2) log(4/3)
    verify(
        SixVertex(; a=1.0, b=1.0, c=1.0),
        ResidualEntropy(),
        Infinite();
        route=:second_closed_form,
        independent=(3 / 2) * log(4 / 3),
        agree_within=1e-9,
        refs=["Lieb 1967: square-ice residual entropy S = (3/2) log(4/3) ≈ 0.4315"],
    )

    # Ferroelectric phase (Δ > 1, e.g. a large): f = -log(max(a,b)), S = 0
    verify(
        SixVertex(; a=3.0, b=1.0, c=1.0),
        ResidualEntropy(),
        Infinite();
        route=:limiting_case,
        independent=0.0,
        agree_within=1e-10,
        refs=["Ferroelectric phase (Δ>1): frozen, residual entropy S = 0"],
    )
end
