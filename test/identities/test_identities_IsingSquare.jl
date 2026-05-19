using Test
using QAtlas

# Self-validation of IsingSquare's finite-T thermodynamics through the
# `verify_thermodynamic_identities` harness.  IsingSquare is a classical
# 2D model — it has no `h` field, so the field-perturbation identity
# (`MAGNETIZATION_X_FROM_FREE_ENERGY`) skips; the three universal
# thermodynamic checks (Gibbs, c_v from ε, c_v from s) all run.

@testset "IsingSquare PBC(4×4) — DEFAULT thermo identities" begin
    # IsingSquare's `Energy{:per_site}` / `SpecificHeat` go through a
    # central-difference stack on `log Z`, so the c_v identities pick up
    # the `O(δ²)` truncation error of the chained derivative; relax the
    # harness atol to `1e-4` to accept those.  Gibbs (algebraic) still
    # holds at machine precision.
    #
    # `IsingSquare(; Lx, Ly)` carries the lattice extents on the struct
    # so the harness `fetch(model, ..., bc; beta=…)` (no Lx/Ly kwargs)
    # picks them up via the `Lx=m.Lx, Ly=m.Ly` defaults.
    model = IsingSquare(; J=1.0, Lx=4, Ly=4)
    βs = [0.1, 0.3, 0.6]
    results = verify_thermodynamic_identities(model, PBC(0); βs=βs, atol=1e-4)

    # 4 default identities × 3 βs.  m_x identity is :skipped (no h
    # field on IsingSquare).
    @test length(results) == 12
    @test all(r.status !== :fail for r in results)
    @test count(r -> r.status === :skipped, results) == 3
    @test any(occursin("Gibbs", r.identity) for r in results)
    @test any(occursin("c_v", r.identity) for r in results)
end

@testset "IsingSquare Infinite — DEFAULT thermo identities (off-T_c)" begin
    # T_c ≈ β_c = log(1+√2)/2 ≈ 0.4407.  Stay clear of it: at the
    # critical point Onsager c_v diverges, and the Gibbs identity is
    # dominated by the c_v residual which blows up under any finite δ
    # central difference.
    model = IsingSquare(; J=1.0)
    βs = [0.15, 0.30, 0.65]
    results = verify_thermodynamic_identities(model, Infinite(); βs=βs, atol=1e-3)

    @test length(results) == 12
    @test all(r.status !== :fail for r in results)
    @test count(r -> r.status === :skipped, results) == 3
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "IsingSquare identities — verification cards" begin
    # Onsager 1944 critical temperature: sinh(2 βc J) = 1.
    for J in (1.0, 2.0)
        verify(
            IsingSquare(; J=J),
            CriticalTemperature(),
            Infinite();
            route=:second_closed_form,
            independent=2 * J / log(1 + sqrt(2)),
            agree_within=1e-10,
            refs=["Onsager 1944: Tc = 2J / log(1+√2)"],
        )
    end
end
