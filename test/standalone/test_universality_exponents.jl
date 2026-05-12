# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: universality class critical exponents
#
# Verify exact (rational) and numerical exponents, and check standard
# scaling relations.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

# ─────────────── Helper: check scaling relations on exact exponents ───────────

function check_scaling_relations(e; d::Int=0)
    @test e.α + 2 * e.β + e.γ == 2                 # Rushbrooke
    @test e.γ == e.β * (e.δ - 1)                    # Widom
    @test e.γ == e.ν * (2 - e.η)                    # Fisher
    if d > 0
        @test 2 - e.α ≈ d * e.ν atol = 1e-10       # Josephson (hyperscaling)
    end
end

# `atol = 0.01` is a *physics-imposed* floor, not laziness — see #118 audit.
#
# The 3D bootstrap exponents we ship (Ising / XY / Heisenberg) inherit the
# precision of their source data (KPSDV 2016 etc.).  Worst-case residuals
# of `α + 2β + γ - 2`, `γ - β(δ-1)`, `γ - ν(2-η)`, `2 - α - dν` measured
# at PR time:
#
#   3D Ising       — max ≈ 1.0e-5  (could tighten to atol = 1e-4)
#   3D XY          — max ≈ 7.1e-5  (could tighten to atol = 1e-3)
#   3D Heisenberg  — max ≈ 4.5e-4  (cannot tighten below atol ≈ 1e-3)
#
# Heisenberg sets the floor: tightening per call site would only shave a
# factor or two off the loosest class while opening every numerical-data
# refresh to spurious failures.  The single `atol = 0.01` here gives ~20×
# margin uniformly and stays robust to upstream value updates.  This is
# the (a)/(b)/(c)-style classification from #118: a deliberately soft
# tolerance with a one-line justification.
function check_scaling_relations_approx(e; d::Int=0)
    @test e.α + 2 * e.β + e.γ ≈ 2 atol = 0.01      # Rushbrooke
    @test e.γ ≈ e.β * (e.δ - 1) atol = 0.01         # Widom
    @test e.γ ≈ e.ν * (2 - e.η) atol = 0.01         # Fisher
    if d > 0
        @test 2 - e.α ≈ d * e.ν atol = 0.01         # Josephson
    end
end

# ═══════════════════ Exact universality classes (Rational) ═══════════════════

@testset "Universality: 2D Ising (exact)" begin
    e = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2)
    @test e.β == 1 // 8
    @test e.ν == 1 // 1
    @test e.γ == 7 // 4
    @test e.η == 1 // 4
    @test e.δ == 15 // 1
    @test e.α == 0 // 1
    @test e.c == 1 // 2
    check_scaling_relations(e; d=2)
end

@testset "Universality: backward compat Ising2D()" begin
    e_old = QAtlas.fetch(Ising2D(), CriticalExponents())
    e_new = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2)
    @test e_old == e_new
end

@testset "Universality: Mean-field (exact)" begin
    e = QAtlas.fetch(MeanField(), CriticalExponents())
    @test e.β == 1 // 2
    @test e.ν == 1 // 2
    @test e.γ == 1 // 1
    @test e.η == 0 // 1
    @test e.δ == 3 // 1
    @test e.α == 0 // 1
    check_scaling_relations(e)
end

@testset "Universality: Ising d≥4 = Mean-field" begin
    e_mf = QAtlas.fetch(MeanField(), CriticalExponents())
    for d in [4, 5, 100]
        e = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=d)
        @test e == e_mf
    end
end

@testset "Universality: 2D Percolation (exact)" begin
    e = QAtlas.fetch(Universality(:Percolation), CriticalExponents(); d=2)
    @test e.α == -2 // 3
    @test e.β == 5 // 36
    @test e.γ == 43 // 18
    @test e.δ == 91 // 5
    @test e.ν == 4 // 3
    @test e.η == 5 // 24
    check_scaling_relations(e; d=2)
end

@testset "Universality: 3-state Potts d=2 (exact)" begin
    e = QAtlas.fetch(Universality(:Potts3), CriticalExponents(); d=2)
    @test e.β == 1 // 9
    @test e.ν == 5 // 6
    @test e.η == 4 // 15
    @test e.δ == 14 // 1
    check_scaling_relations(e; d=2)
end

@testset "Universality: 4-state Potts d=2 (exact)" begin
    e = QAtlas.fetch(Universality(:Potts4), CriticalExponents(); d=2)
    @test e.β == 1 // 12
    @test e.ν == 2 // 3
    @test e.η == 1 // 4
    @test e.δ == 15 // 1
    check_scaling_relations(e; d=2)
end

@testset "Universality: KPZ 1+1D (exact growth exponents)" begin
    e = QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=1)
    @test e.β_growth == 1 // 3
    @test e.α_rough == 1 // 2
    @test e.z == 3 // 2
    @test e.α_rough + e.z == 2              # Galilean invariance
    @test e.β_growth == e.α_rough / e.z     # β = α / z
end

@testset "Universality: KPZ1D() backward compat" begin
    e_old = QAtlas.fetch(KPZ1D(), CriticalExponents())
    e_new = QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=1)
    @test e_old == e_new
end

# ═══════════════ Numerical universality classes (Float64 + _err) ═════════════

@testset "Universality: 3D Ising (conformal bootstrap)" begin
    e = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=3)
    @test 0.10 < e.β < 0.35
    @test 0.60 < e.ν < 0.65
    @test e.β_err > 0
    @test e.ν_err > 0
    check_scaling_relations_approx(e; d=3)
end

@testset "Universality: 3D XY (conformal bootstrap)" begin
    e = QAtlas.fetch(Universality(:XY), CriticalExponents(); d=3)
    @test 0.34 < e.β < 0.36
    @test 0.66 < e.ν < 0.68
    @test e.β_err > 0
    check_scaling_relations_approx(e; d=3)
end

@testset "Universality: 3D Heisenberg" begin
    e = QAtlas.fetch(Universality(:Heisenberg), CriticalExponents(); d=3)
    @test 0.36 < e.β < 0.38
    @test 0.70 < e.ν < 0.72
    @test e.β_err > 0
    check_scaling_relations_approx(e; d=3)
end

@testset "Universality: XY d=2 is BKT" begin
    e = QAtlas.fetch(Universality(:XY), CriticalExponents(); d=2)
    @test e.η == 1 // 4
end

@testset "Universality: Heisenberg d=2 → Mermin-Wagner error" begin
    @test_throws ErrorException QAtlas.fetch(
        Universality(:Heisenberg), CriticalExponents(); d=2
    )
end

@testset "Universality: 3D Percolation (numerical)" begin
    e = QAtlas.fetch(Universality(:Percolation), CriticalExponents(); d=3)
    @test 0.40 < e.β < 0.45
    @test 0.85 < e.ν < 0.90
    @test e.β_err > 0
end

# ═══════════════ KPZ higher dimensions (numerical) ════════════════════════════

@testset "Universality: KPZ 2+1D (Pagnani–Parisi 2015 numerical)" begin
    e = QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=2)
    # Pagnani & Parisi, PRE 92, 010101(R) (2015)
    @test 0.235 < e.β_growth < 0.245
    @test 0.385 < e.α_rough < 0.400
    @test 1.600 < e.z < 1.625
    @test e.β_growth_err > 0
    @test e.α_rough_err > 0
    @test e.z_err > 0
    # Galilean invariance α + z = 2 should hold within combined error.
    # Use a generous 3 σ to absorb correlated estimation bias.
    @test abs(e.α_rough + e.z - 2) < 3 * (e.α_rough_err + e.z_err)
end

@testset "Universality: KPZ 3+1D (Kelling–Ódor 2011 numerical)" begin
    e = QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=3)
    # Kelling & Ódor, PRE 84, 061150 (2011)
    @test 0.16 < e.β_growth < 0.20
    @test 0.29 < e.α_rough < 0.33
    @test 1.49 < e.z < 1.53
    @test e.β_growth_err > 0
    # Note: Galilean invariance is NOT satisfied within quoted errors at d=3
    # (α + z ≈ 1.82 vs 2.0).  This reflects open issues in the d≥3 KPZ
    # literature; the stored values are best-numerical estimates, not a
    # rigorous reference.  Therefore: no α + z = 2 assertion here.
end

@testset "Universality: KPZ d≥4 errors" begin
    @test_throws ErrorException QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=4)
    @test_throws ErrorException QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=5)
end
