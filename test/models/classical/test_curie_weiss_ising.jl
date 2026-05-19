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
