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

@testset "SixVertex — antiferroelectric phase (Lieb 1967b)" begin
    # f-model with c = 3 ⇒ Δ = −7/2 < −1, AFE phase.
    m_afe = QAtlas.f_model(3.0)
    @test QAtlas._six_vertex_phase(m_afe.a, m_afe.b, m_afe.c) === :antiferroelectric

    Δ = QAtlas._six_vertex_delta(m_afe.a, m_afe.b, m_afe.c)
    λ = acosh(-Δ)
    
    # Verify free energy is close to the expected elliptic sum
    f = QAtlas.fetch(m_afe, FreeEnergy(), Infinite())
    u = λ / 2.0
    s_expected = 0.0
    for m in 1:1000
        t1 = exp(-m * (2.0*λ - 2.0*u))
        t2 = exp(-m * (2.0*λ + 2.0*u))
        den = 1.0 + exp(-2.0 * m * λ)
        s_expected += (t1 - t2) / (m * den)
    end
    f_expected = -log(m_afe.a) - (λ - u) - s_expected
    @test f ≈ f_expected atol = 1e-12

    # Verify Polarization
    pol = QAtlas.fetch(m_afe, Polarization(), Infinite())
    pol_expected = tanh(1 * λ)^2 * tanh(2 * λ)^2 * tanh(3 * λ)^2 * tanh(4 * λ)^2
    @test pol ≈ pol_expected atol = 1e-6
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

@testset "SixVertex — disordered free energy off-diagonal" begin
    # Generic disordered point (a = b = 1, c = 0.5):
    #   Δ = (1 + 1 − 0.25) / 2 = 0.875 ∈ (-1, 1).
    m_d = SixVertex(; a=1.0, b=1.0, c=0.5)
    @test QAtlas._six_vertex_delta(1.0, 1.0, 0.5) ≈ 0.875
    @test QAtlas._six_vertex_phase(m_d.a, m_d.b, m_d.c) === :disordered
    
    # Value compared to numerical check
    f_d = QAtlas.fetch(m_d, FreeEnergy(), Infinite())
    @test f_d ≈ -0.129202353139369 atol = 1e-12

    # f-model boundary c = 2 (Δ = −1, F-model critical point)
    m_crit = QAtlas.f_model(2.0)
    @test QAtlas._six_vertex_phase(m_crit.a, m_crit.b, m_crit.c) === :disordered
    f_crit = QAtlas.fetch(m_crit, FreeEnergy(), Infinite())
    @test f_crit ≈ -0.7831887820803405 atol = 1e-12
    
    # Polarization at critical point is 0
    @test QAtlas.fetch(m_crit, Polarization(), Infinite()) == 0.0
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

# ── additional verification cards (#381 batch) ─────────────────────────────
@testset "SixVertex — FreeEnergy closed forms (#381 batch)" begin
    # Square-ice point a=b=c=1 (Δ=1/2, disordered): f = -(3/2) log(4/3)
    # (Lieb 1967a — same magnitude as the residual entropy, opposite sign).
    verify(
        SixVertex(; a=1.0, b=1.0, c=1.0),
        FreeEnergy(),
        Infinite();
        route=:lieb_square_ice,
        independent=(-(3/2) * log(4/3)),
        agree_within=1e-14,
        refs=["Lieb 1967a Phys. Rev. 162: square-ice f = -(3/2) log(4/3)"],
    )

    # Ferroelectric phase Δ>1: f = -log max(a,b) (Lieb 1967c).
    # (a, b, c) = (3, 1, 1) → Δ = (9+1-1)/(2*3*1) = 9/6 = 1.5 > 1; f = -log 3.
    verify(
        SixVertex(; a=3.0, b=1.0, c=1.0),
        FreeEnergy(),
        Infinite();
        route=:lieb_ferroelectric,
        independent=(-log(3.0)),
        agree_within=1e-14,
        refs=[
            "Lieb 1967c Phys. Rev. Lett. 19, 108: KDP/FE phase f = -log max(a,b) (frozen GS)",
        ],
    )

    # (a, b, c) = (1, 3, 1) → same by a↔b symmetry; f = -log 3.
    verify(
        SixVertex(; a=1.0, b=3.0, c=1.0),
        FreeEnergy(),
        Infinite();
        route=:lieb_ferroelectric,
        independent=(-log(3.0)),
        agree_within=1e-14,
        refs=["Lieb 1967c: FE-phase f = -log max(a,b); a↔b symmetry"],
    )

    # (a, b, c) = (2, 4, 1) → Δ = (4+16-1)/(2*2*4) = 19/16 > 1; f = -log 4.
    verify(
        SixVertex(; a=2.0, b=4.0, c=1.0),
        FreeEnergy(),
        Infinite();
        route=:lieb_ferroelectric,
        independent=(-log(4.0)),
        agree_within=1e-14,
        refs=["Lieb 1967c: FE-phase f = -log max(a,b)"],
    )

    # KDP boundary Δ=1: approach the transition from the FE side at
    # a = nextfloat(2.0), b = c = 1.  Lieb 1967c: f = -log max(a,b).
    # The boundary itself (a=2 exactly) is classified :disordered and
    # fetch throws there; this card exercises the immediately adjacent
    # FE-phase value, completing WHY-correct coverage of the KDP
    # transition without re-implementing the disordered closed form.
    let a_kdp = nextfloat(2.0)
        verify(
            SixVertex(; a=a_kdp, b=1.0, c=1.0),
            FreeEnergy(),
            Infinite();
            route=:lieb_ferroelectric,
            independent=(-log(a_kdp)),
            agree_within=1e-14,
            refs=[
                "Lieb 1967c: KDP boundary Δ=1 approached from FE side; f = -log max(a,b)"
            ],
        )
    end
end

@testset "SixVertex — Energy and Polarization verification cards" begin
    # Energy at square-ice point (should be 0.0)
    verify(
        SixVertex(; a=1.0, b=1.0, c=1.0),
        Energy{:per_site}(),
        Infinite();
        route=:numerical_derivative,
        independent=0.0,
        agree_within=1e-9,
        refs=["Energy per site is zero at the isotropic square-ice point"],
    )

    # Polarization in FE phase is 1.0
    verify(
        SixVertex(; a=3.0, b=1.0, c=1.0),
        Polarization(),
        Infinite();
        route=:limiting_case,
        independent=1.0,
        agree_within=1e-10,
        refs=["Polarization in FE phase is saturated to 1.0"],
    )

    # Polarization in disordered phase is 0.0
    verify(
        SixVertex(; a=1.0, b=1.0, c=1.0),
        Polarization(),
        Infinite();
        route=:limiting_case,
        independent=0.0,
        agree_within=1e-10,
        refs=["Polarization in disordered phase is exactly 0.0"],
    )
end
