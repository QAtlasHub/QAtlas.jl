# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: universality class critical exponents
#
# Migrated from pure-legacy @test to verify()-first (PR #449 phase B,
# zero-legacy end-state). Exact rational exponents become verify() cards
# via subject_extract = e -> getproperty(e, field); numerical bounded
# ranges, _err positivity, scaling-relation helper identities, KPZ
# higher-d numerical estimates, and Mermin-Wagner / KPZ d>=4 error paths
# stay raw @test (verify() is scalar-only and cannot represent
# inequalities or error-throwing dispatches).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

# ─────────────── Helper: check scaling relations on exact exponents ───────────

function check_scaling_relations(e; d::Int=0)
    @test e.α + 2 * e.β + e.γ == 2                 # Rushbrooke
    @test e.γ == e.β * (e.δ - 1)                    # Widom
    @test e.γ == e.ν * (2 - e.η)                    # Fisher
    if d > 0
        @test 2 - e.α ≈ d * e.ν atol = 1e-10        # Josephson
    end
end

function check_scaling_relations_approx(e; d::Int=0)
    @test e.α + 2 * e.β + e.γ ≈ 2 atol = 0.01
    @test e.γ ≈ e.β * (e.δ - 1) atol = 0.01
    @test e.γ ≈ e.ν * (2 - e.η) atol = 0.01
    if d > 0
        @test 2 - e.α ≈ d * e.ν atol = 0.01
    end
end

# ═══════════════════ Exact universality classes (verify) ═════════════════════

@testset "Universality: 2D Ising (exact)" begin
    e = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2)
    for (field, lit) in (
        (:β, 1 // 8),
        (:ν, 1 // 1),
        (:γ, 7 // 4),
        (:η, 1 // 4),
        (:δ, 15 // 1),
        (:α, 0 // 1),
        (:c, 1 // 2),
    )
        verify(
            Universality(:Ising),
            CriticalExponents(),
            Infinite();
            route=:literature_value,
            independent=lit,
            agree_within=0,
            at=["d=2", "field=$(field)"],
            refs=["Onsager 1944 / Pfeuty 1970: 2D Ising exact exponent $(field) = $(lit)"],
            fetch_kw=(; d=2),
            subject_extract=e -> getproperty(e, field),
        )
    end
    check_scaling_relations(e; d=2)
end

@testset "Universality: backward compat Ising2D()" begin
    e_old = QAtlas.fetch(Ising2D(), CriticalExponents())
    e_new = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2)
    @test e_old == e_new
end

@testset "Universality: Mean-field (exact)" begin
    for (field, lit) in
        ((:β, 1 // 2), (:ν, 1 // 2), (:γ, 1 // 1), (:η, 0 // 1), (:δ, 3 // 1), (:α, 0 // 1))
        verify(
            MeanField(),
            CriticalExponents(),
            Infinite();
            route=:literature_value,
            independent=lit,
            agree_within=0,
            at=["field=$(field)"],
            refs=["Landau mean-field / upper critical d: $(field) = $(lit)"],
            subject_extract=e -> getproperty(e, field),
        )
    end
    check_scaling_relations(QAtlas.fetch(MeanField(), CriticalExponents()))
end

@testset "Universality: Ising d≥4 = Mean-field" begin
    e_mf = QAtlas.fetch(MeanField(), CriticalExponents())
    for d in (4, 5, 100)
        e = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=d)
        @test e == e_mf
    end
end

@testset "Universality: 2D Percolation (exact)" begin
    for (field, lit) in (
        (:α, -2 // 3),
        (:β, 5 // 36),
        (:γ, 43 // 18),
        (:δ, 91 // 5),
        (:ν, 4 // 3),
        (:η, 5 // 24),
    )
        verify(
            Universality(:Percolation),
            CriticalExponents(),
            Infinite();
            route=:literature_value,
            independent=lit,
            agree_within=0,
            at=["d=2", "field=$(field)"],
            refs=[
                "Stauffer-Aharony / den Nijs 1979: 2D Percolation exact $(field) = $(lit)"
            ],
            fetch_kw=(; d=2),
            subject_extract=e -> getproperty(e, field),
        )
    end
    e = QAtlas.fetch(Universality(:Percolation), CriticalExponents(); d=2)
    check_scaling_relations(e; d=2)
end

@testset "Universality: 3-state Potts d=2 (exact)" begin
    for (field, lit) in ((:β, 1 // 9), (:ν, 5 // 6), (:η, 4 // 15), (:δ, 14 // 1))
        verify(
            Universality(:Potts3),
            CriticalExponents(),
            Infinite();
            route=:literature_value,
            independent=lit,
            agree_within=0,
            at=["d=2", "field=$(field)"],
            refs=["Dotsenko 1984 / DFMS §7.4: 2D 3-state Potts $(field) = $(lit)"],
            fetch_kw=(; d=2),
            subject_extract=e -> getproperty(e, field),
        )
    end
    e = QAtlas.fetch(Universality(:Potts3), CriticalExponents(); d=2)
    check_scaling_relations(e; d=2)
end

@testset "Universality: 4-state Potts d=2 (exact)" begin
    for (field, lit) in ((:β, 1 // 12), (:ν, 2 // 3), (:η, 1 // 4), (:δ, 15 // 1))
        verify(
            Universality(:Potts4),
            CriticalExponents(),
            Infinite();
            route=:literature_value,
            independent=lit,
            agree_within=0,
            at=["d=2", "field=$(field)"],
            refs=[
                "DFMS §12.3: 2D 4-state Potts (marginal compact boson) $(field) = $(lit)"
            ],
            fetch_kw=(; d=2),
            subject_extract=e -> getproperty(e, field),
        )
    end
    e = QAtlas.fetch(Universality(:Potts4), CriticalExponents(); d=2)
    check_scaling_relations(e; d=2)
end

@testset "Universality: KPZ 1+1D (exact growth exponents)" begin
    for (field, lit) in ((:β_growth, 1 // 3), (:α_rough, 1 // 2), (:z, 3 // 2))
        verify(
            Universality(:KPZ),
            GrowthExponents(),
            Infinite();
            route=:literature_value,
            independent=lit,
            agree_within=0,
            at=["d=1", "field=$(field)"],
            refs=["KPZ 1+1D exact (Kardar-Parisi-Zhang 1986): $(field) = $(lit)"],
            fetch_kw=(; d=1),
            subject_extract=e -> getproperty(e, field),
        )
    end
    e = QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=1)
    @test e.α_rough + e.z == 2              # Galilean invariance (derived identity, raw)
    @test e.β_growth == e.α_rough / e.z     # β = α / z (derived, raw)
end

@testset "Universality: KPZ1D() backward compat" begin
    e_old = QAtlas.fetch(KPZ1D(), CriticalExponents())
    e_new = QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=1)
    @test e_old == e_new
end

# ═══════════════ Numerical universality classes — kept raw ════════════════════

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
    verify(
        Universality(:XY),
        CriticalExponents(),
        Infinite();
        route=:literature_value,
        independent=1 // 4,
        agree_within=0,
        at=["d=2", "field=η"],
        refs=["Kosterlitz 1974 BKT: η = 1/4 at the universal jump"],
        fetch_kw=(; d=2),
        subject_extract=e -> e.η,
    )
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

# ═══════════════ KPZ higher dimensions (numerical) ────────────────────────────

@testset "Universality: KPZ 2+1D (Pagnani-Parisi 2015 numerical)" begin
    e = QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=2)
    @test 0.235 < e.β_growth < 0.245
    @test 0.385 < e.α_rough < 0.400
    @test 1.600 < e.z < 1.625
    @test e.β_growth_err > 0
    @test e.α_rough_err > 0
    @test e.z_err > 0
    @test abs(e.α_rough + e.z - 2) < 3 * (e.α_rough_err + e.z_err)
end

@testset "Universality: KPZ 3+1D (Kelling-Ódor 2011 numerical)" begin
    e = QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=3)
    @test 0.16 < e.β_growth < 0.20
    @test 0.29 < e.α_rough < 0.33
    @test 1.49 < e.z < 1.53
    @test e.β_growth_err > 0
end

@testset "Universality: KPZ d≥4 errors" begin
    @test_throws ErrorException QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=4)
    @test_throws ErrorException QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=5)
end
