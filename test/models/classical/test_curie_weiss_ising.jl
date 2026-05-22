# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: CurieWeissIsing — complete-graph mean-field Ising.
#
# Verifies:
#   * T_c = J (J > 0); T_c = 0 (J ≤ 0)
#   * Spontaneous magnetisation m*(β) solves m = tanh(βJ m) with the
#     positive root being correct to 1e-12; m = 0 in the paramagnetic
#     phase β J ≤ 1
#   * Landau small-(T_c - T) expansion m*² ≈ 3 (1 - T/T_c) at leading
#     order — independent check that the Newton solver found the
#     physical branch
#   * Mean-field exponent β = 1/2: log m* vs log(1 - T/T_c) has slope
#     1/2 close to T_c
#   * β ≤ 0 throws DomainError
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "CurieWeissIsing — CriticalTemperature" begin
    @test QAtlas.fetch(CurieWeissIsing(; J=1.0), CriticalTemperature(), Infinite()) ≈ 1.0
    @test QAtlas.fetch(CurieWeissIsing(; J=2.5), CriticalTemperature(), Infinite()) ≈ 2.5
    # J ≤ 0: no FM order, T_c = 0
    @test QAtlas.fetch(CurieWeissIsing(; J=0.0), CriticalTemperature(), Infinite()) == 0.0
    @test QAtlas.fetch(CurieWeissIsing(; J=-1.5), CriticalTemperature(), Infinite()) == 0.0
end

@testset "CurieWeissIsing — SpontaneousMagnetization paramagnetic phase" begin
    m = CurieWeissIsing(; J=1.0)
    # T > T_c (βJ < 1): m* = 0
    for β in (0.1, 0.5, 0.99, 1.0)
        @test QAtlas.fetch(m, SpontaneousMagnetization(), Infinite(); beta=β) == 0.0
    end
end

@testset "CurieWeissIsing — SpontaneousMagnetization self-consistency" begin
    # Verify the returned m* satisfies m = tanh(βJ m) to high precision.
    m = CurieWeissIsing(; J=1.0)
    for β in (1.01, 1.1, 1.5, 2.0, 3.0, 5.0)
        m_star = QAtlas.fetch(m, SpontaneousMagnetization(), Infinite(); beta=β)
        @test m_star > 0
        @test tanh(β * m_star) ≈ m_star atol = 1e-12
    end
end

@testset "CurieWeissIsing — Landau β = 1/2 critical exponent" begin
    # Just below T_c with reduced temperature t = (T_c - T)/T_c, the
    # Landau theory predicts m*(t) ≈ √(3 t) at leading order.
    m = CurieWeissIsing(; J=1.0)
    ts = [1e-4, 5e-5]                          # T_c - T  (with T_c = 1)
    βs = 1 ./ (1 .- ts)                         # β = 1/T = 1/(1 - t)
    mstars = [QAtlas.fetch(m, SpontaneousMagnetization(), Infinite(); beta=β) for β in βs]
    @test mstars[1] ≈ sqrt(3 * ts[1]) rtol = 5e-3
    @test mstars[2] ≈ sqrt(3 * ts[2]) rtol = 5e-3
    # Slope in log-log space ≈ β_exp = 1/2.
    slope = log(mstars[1] / mstars[2]) / log(ts[1] / ts[2])
    @test slope ≈ 1 / 2 atol = 5e-3
end

@testset "CurieWeissIsing — SpontaneousMagnetization deep ordered phase" begin
    # As T → 0, m*(β) → 1.  Pin a value at β = 10 J (well below T_c).
    m = CurieWeissIsing(; J=1.0)
    m_star = QAtlas.fetch(m, SpontaneousMagnetization(), Infinite(); beta=10.0)
    @test 0.99 < m_star < 1.0
    @test tanh(10 * m_star) ≈ m_star atol = 1e-12
end

@testset "CurieWeissIsing — m* depends on β and J only through βJ" begin
    m1 = CurieWeissIsing(; J=1.0)
    m2 = CurieWeissIsing(; J=2.0)
    val1 = QAtlas.fetch(m1, SpontaneousMagnetization(), Infinite(); beta=1.5)
    val2 = QAtlas.fetch(m2, SpontaneousMagnetization(), Infinite(); beta=0.75)
    @test val1 ≈ val2 atol = 1e-12
end

@testset "CurieWeissIsing — DomainError on non-positive beta" begin
    m = CurieWeissIsing(; J=1.0)
    @test_throws DomainError QAtlas.fetch(m, SpontaneousMagnetization(), Infinite(); beta=0)
    @test_throws DomainError QAtlas.fetch(
        m, SpontaneousMagnetization(), Infinite(); beta=-1.0
    )
end

@testset "CurieWeissIsing — J ≤ 0 paramagnet at all β" begin
    for J in (0.0, -1.0)
        m = CurieWeissIsing(; J=J)
        @test QAtlas.fetch(m, SpontaneousMagnetization(), Infinite(); beta=10.0) == 0.0
    end
end

@testset "CurieWeissIsing — CriticalExponents delegate to MeanField (Phase 2)" begin
    m = CurieWeissIsing()
    exp = QAtlas.fetch(m, CriticalExponents(), Infinite())
    @test exp.α == 0
    @test exp.β == 1 // 2
    @test exp.γ == 1
    @test exp.δ == 3
    @test exp.ν == 1 // 2
    @test exp.η == 0
    # Delegation invariant
    @test exp == QAtlas.fetch(QAtlas.MeanField(), CriticalExponents())
    # Hyperscaling: Rushbrooke α + 2β + γ = 2; Widom γ = β(δ − 1)
    @test exp.α + 2 * exp.β + exp.γ == 2
    @test exp.γ == exp.β * (exp.δ - 1)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "CurieWeissIsing — verification cards" begin
    # Mean-field critical temperature Tc = J (independent: the
    # self-consistency m = tanh(βJm) linearises to βc J = 1 => Tc = J).
    for J in (0.5, 1.0, 2.5)
        verify(
            CurieWeissIsing(; J=J),
            CriticalTemperature(),
            Infinite();
            route=:second_closed_form,
            independent=J,
            agree_within=1e-12,
            refs=["Mean-field: linearised self-consistency gives βc J = 1 => Tc = J"],
        )
    end

    # Spontaneous magnetization solves m = tanh(βJm); independently
    # re-solved here by fixed-point iteration (not from src).
    let J = 1.0, β = 2.0
        m = 0.9
        for _ in 1:200
            m = tanh(β * J * m)
        end
        verify(
            CurieWeissIsing(; J=J),
            SpontaneousMagnetization(),
            Infinite();
            route=:second_closed_form,
            fetch_kw=(; beta=β),
            independent=m,
            agree_within=1e-8,
            refs=["Curie-Weiss self-consistency m = tanh(βJm), independent fixed point"],
        )
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Mean-field thermodynamics — verify cards covering both h = 0 and h ≠ 0
# ═══════════════════════════════════════════════════════════════════════════════

# ── h = 0 cards: cross-check derived quantities against independent limits
#    that don't re-run the SCE solver ──────────────────────────────────────
@testset "CurieWeissIsing — FreeEnergy at h=0, T > T_c: paramagnet -log(2)/β" begin
    for (J, β) in ((1.0, 0.5), (1.0, 0.9), (2.0, 0.3))   # βJ < 1
        verify(
            CurieWeissIsing(; J=J, h=0.0),
            FreeEnergy(),
            Infinite();
            route=:limiting_case,
            fetch_kw=(; beta=β),
            independent=-log(2) / β,
            agree_within=1e-12,
            refs=["Landau-Lifshitz §149: paramagnet (m*=0) ⇒ f = -log(2)/β per site"],
        )
    end
end

@testset "CurieWeissIsing — Energy/per_site at h=0, T > T_c: u = 0" begin
    for (J, β) in ((1.0, 0.5), (1.0, 0.9), (2.0, 0.3))
        verify(
            CurieWeissIsing(; J=J, h=0.0),
            Energy{:per_site}(),
            Infinite();
            route=:limiting_case,
            fetch_kw=(; beta=β),
            independent=0.0,
            agree_within=1e-14,
            refs=["Landau-Lifshitz §149: paramagnet ⇒ u = -Jm*²/2 = 0"],
        )
    end
end

@testset "CurieWeissIsing — ThermalEntropy at h=0, T > T_c: s = log 2" begin
    for (J, β) in ((1.0, 0.5), (1.0, 0.9), (2.0, 0.3))
        verify(
            CurieWeissIsing(; J=J, h=0.0),
            ThermalEntropy(),
            Infinite();
            route=:limiting_case,
            fetch_kw=(; beta=β),
            independent=log(2.0),
            agree_within=1e-14,
            refs=["Landau-Lifshitz §149: paramagnet ⇒ s = log 2 per site"],
        )
    end
end

@testset "CurieWeissIsing — SpecificHeat at h=0, T > T_c: c_v = 0" begin
    for (J, β) in ((1.0, 0.5), (1.0, 0.9), (2.0, 0.3))
        verify(
            CurieWeissIsing(; J=J, h=0.0),
            SpecificHeat(),
            Infinite();
            route=:limiting_case,
            fetch_kw=(; beta=β),
            independent=0.0,
            agree_within=1e-14,
            refs=["Landau-Lifshitz §149: paramagnet ⇒ c_v = 0 (Landau jump at T_c⁻)"],
        )
    end
end

@testset "CurieWeissIsing — SusceptibilityZZ at h=0, T > T_c: Curie-Weiss law" begin
    # χ(β,J,h=0) = β / (1 - βJ) for βJ < 1
    for (J, β) in ((1.0, 0.5), (1.0, 0.9), (2.0, 0.3))
        verify(
            CurieWeissIsing(; J=J, h=0.0),
            SusceptibilityZZ(),
            Infinite();
            route=:second_closed_form,
            fetch_kw=(; beta=β),
            independent=β / (1 - β * J),
            agree_within=1e-12,
            refs=["Curie-Weiss law (Landau-Lifshitz §149): χ = β/(1-βJ) above T_c"],
        )
    end
end

# ── h ≠ 0 cards: J = 0 reduces to a single spin in a transverse field,
#    where every quantity has a textbook closed form independent of any
#    SCE. ────────────────────────────────────────────────────────────────
@testset "CurieWeissIsing — J=0, h≠0: single spin in field, all 5 quantities" begin
    # Reference: single spin in field h with H = -h σ at inverse temp β.
    #   u = -h tanh(βh)
    #   f = -log(2 cosh βh) / β
    #   s = log(2 cosh βh) - βh tanh βh
    #   c_v = (βh sech βh)²
    #   χ = β sech²(βh)
    # Independent of the CurieWeiss SCE since J=0 short-circuits the
    # mean-field branch.
    for h in (0.1, 0.5, 1.5)
        for β in (0.3, 1.0, 3.0)
            m = CurieWeissIsing(; J=0.0, h=h)
            verify(
                m,
                FreeEnergy(),
                Infinite();
                route=:second_closed_form,
                fetch_kw=(; beta=β),
                independent=-log(2 * cosh(β * h)) / β,
                agree_within=1e-12,
                refs=["Single spin in field h: f = -β⁻¹ log(2 cosh βh)"],
            )
            verify(
                m,
                Energy{:per_site}(),
                Infinite();
                route=:second_closed_form,
                fetch_kw=(; beta=β),
                independent=-h * tanh(β * h),
                agree_within=1e-12,
                refs=["Single spin in field h: u = -h tanh(βh)"],
            )
            verify(
                m,
                ThermalEntropy(),
                Infinite();
                route=:second_closed_form,
                fetch_kw=(; beta=β),
                independent=log(2 * cosh(β * h)) - β * h * tanh(β * h),
                agree_within=1e-12,
                refs=["Single spin in field h: s = log(2 cosh βh) - βh tanh(βh)"],
            )
            verify(
                m,
                SpecificHeat(),
                Infinite();
                route=:second_closed_form,
                fetch_kw=(; beta=β),
                independent=(β * h * sech(β * h))^2,
                agree_within=1e-12,
                refs=["Single spin in field h: c_v = (βh sech βh)²"],
            )
            verify(
                m,
                SusceptibilityZZ(),
                Infinite();
                route=:second_closed_form,
                fetch_kw=(; beta=β),
                independent=β * sech(β * h)^2,
                agree_within=1e-12,
                refs=["Single spin in field h: χ = β sech²(βh)"],
            )
        end
    end
end

# ── h ≠ 0 cards: T → 0 saturation with J > 0 ─────────────────────────────
@testset "CurieWeissIsing — T → 0 saturation with h > 0: u → -J/2 - h" begin
    # At β → ∞ with J > 0 and h > 0, the SCE forces m* → 1 (saturation
    # by both field and exchange), independent of the precise solver:
    #   u → -J·1²/2 - h·1 = -J/2 - h
    # Verified directly (no need to re-run the bisection).
    for (J, h) in ((1.0, 0.3), (2.0, 0.5), (0.7, 0.1))
        verify(
            CurieWeissIsing(; J=J, h=h),
            Energy{:per_site}(),
            Infinite();
            route=:limiting_case,
            fetch_kw=(; beta=100.0),
            independent=-J / 2 - h,
            agree_within=1e-10,
            refs=[
                "T → 0 saturation: m → 1 (h > 0 + J > 0 selects full alignment) ⇒ u = -J/2 - h",
            ],
        )
    end
end

@testset "CurieWeissIsing — T → 0 saturation with h > 0: s → 0" begin
    for (J, h) in ((1.0, 0.3), (2.0, 0.5))
        verify(
            CurieWeissIsing(; J=J, h=h),
            ThermalEntropy(),
            Infinite();
            route=:limiting_case,
            fetch_kw=(; beta=200.0),
            independent=0.0,
            agree_within=1e-10,     # exact ∞-β limit is 0; residue ~ ULP·β from m = 1 - ε
            refs=["T → 0 saturation: m → 1, fully-ordered ground state ⇒ s = 0"],
        )
    end
end

# ── h ≠ 0 cards: high-T leading-order m ≈ tanh(βh) when βJ ≪ 1 ──────────
@testset "CurieWeissIsing — high-T sanity: u → -h tanh(βh) at βJ ≪ 1" begin
    # At βJ ≪ 1 and h ≠ 0, the Jm* term in the SCE is subleading:
    # m* ≈ tanh(βh) + O(βJ).  So u = -Jm*²/2 - hm* → -h tanh(βh) + O(βJ).
    # We pick β=1e-3, J=0.5, h=0.3 ⇒ βJ=5e-4 (so corrections ~5e-4 relative).
    let β = 1e-3, J = 0.5, h = 0.3
        verify(
            CurieWeissIsing(; J=J, h=h),
            Energy{:per_site}(),
            Infinite();
            route=:limiting_case,
            fetch_kw=(; beta=β),
            independent=-h * tanh(β * h),
            agree_within=1e-3 * abs(h * tanh(β * h)) + 1e-12,
            refs=[
                "βJ → 0 limit: m* → tanh(βh) (J-independent leading order) ⇒ u → -h tanh(βh)",
            ],
        )
    end
end

# ── DomainError and parameter validation ────────────────────────────────
@testset "CurieWeissIsing — DomainError on β ≤ 0 for all thermodynamic quantities" begin
    m = CurieWeissIsing(; J=1.0, h=0.0)
    @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=0.0)
    @test_throws DomainError QAtlas.fetch(m, Energy{:per_site}(), Infinite(); beta=-1.0)
    @test_throws DomainError QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=0.0)
    @test_throws DomainError QAtlas.fetch(m, SpecificHeat(), Infinite(); beta=-2.0)
    @test_throws DomainError QAtlas.fetch(m, SusceptibilityZZ(), Infinite(); beta=0.0)
    @test_throws DomainError QAtlas.fetch(
        m, SpontaneousMagnetization(), Infinite(); beta=0.0
    )
end

@testset "CurieWeissIsing — Gibbs identity cross-check at h ≠ 0" begin
    # s = β(u - f) must hold identically.  Pin a non-trivial point in
    # both ordered (βJ > 1) and disordered (βJ < 1) phase with h ≠ 0.
    for (J, β, h) in ((1.0, 2.0, 0.3), (1.0, 0.5, 0.5), (2.0, 1.0, -0.2))
        m = CurieWeissIsing(; J=J, h=h)
        u = QAtlas.fetch(m, Energy{:per_site}(), Infinite(); beta=β)
        f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
        s = QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=β)
        @test s ≈ β * (u - f) atol = 1e-12
    end
end

# ── Symmetry: u(-h) = u(h), f(-h) = f(h), etc. (all 5 are even in h) ────
@testset "CurieWeissIsing — h → -h symmetry of all quantities" begin
    for (J, β, h) in ((1.0, 2.0, 0.3), (1.0, 0.5, 0.5))
        m_pos = CurieWeissIsing(; J=J, h=h)
        m_neg = CurieWeissIsing(; J=J, h=(-h))
        @test QAtlas.fetch(m_pos, Energy{:per_site}(), Infinite(); beta=β) ≈
            QAtlas.fetch(m_neg, Energy{:per_site}(), Infinite(); beta=β) atol = 1e-12
        @test QAtlas.fetch(m_pos, FreeEnergy(), Infinite(); beta=β) ≈
            QAtlas.fetch(m_neg, FreeEnergy(), Infinite(); beta=β) atol = 1e-12
        @test QAtlas.fetch(m_pos, ThermalEntropy(), Infinite(); beta=β) ≈
            QAtlas.fetch(m_neg, ThermalEntropy(), Infinite(); beta=β) atol = 1e-12
        @test QAtlas.fetch(m_pos, SpecificHeat(), Infinite(); beta=β) ≈
            QAtlas.fetch(m_neg, SpecificHeat(), Infinite(); beta=β) atol = 1e-12
        @test QAtlas.fetch(m_pos, SusceptibilityZZ(), Infinite(); beta=β) ≈
            QAtlas.fetch(m_neg, SusceptibilityZZ(), Infinite(); beta=β) atol = 1e-12
    end
end
