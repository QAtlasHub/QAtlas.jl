# Generated Gibbs identity checks (#698) — the tuple-mode slice of the
# @identity edges (id prefix "identity/gibbs/"), split from the isotropy
# slice so the shard planner can schedule the two independently.  A check
# claimed by NO slice file runs in test_identity_rest.jl — keep the prefix
# lists of the three files in sync.
#
# The hand-written harness (test/util/thermodynamic_identities.jl, exercised
# in test/identities/) stays the reference implementation until the deletion
# criterion in src/identity_registry.jl is met — the anchors below pin the
# generated surface to cover at least the harness's Gibbs hubs.

include("util_run_checks.jl")
using QAtlas: generated_checks

@testset "generated identity checks — :gibbs" begin
    checks = filter(
        c -> startswith(c.id, "identity/gibbs/"), generated_checks(; kinds=(:identity,))
    )
    @test !isempty(checks)

    # Parity anchors: the generated surface covers the hand-written hubs
    # (TFIM at OBC + Infinite) and reaches hubs with NO hand-written identity
    # file — the "new models get identity coverage for free" payoff.
    ids = [c.id for c in checks]
    @test any(startswith("identity/gibbs/TFIM/Infinite"), ids)
    @test any(startswith("identity/gibbs/TFIM/OBC"), ids)
    @test any(startswith("identity/gibbs/XXZ1D/"), ids)

    run_generated_suite(checks; label="generated gibbs checks")
end

# The edge no longer restates `f = ε − T·s`; it asks AbstractQAtlas's
# `FreeEnergyLegendre` to derive `f` (#734 Phase B).  Pin that delegation against
# the closed form it replaced, so an upstream convention change — sign, per-site
# granularity, β-vs-T — fails HERE with a two-line diff instead of silently
# redefining what ":gibbs" asserts for every hub in the atlas.
@testset ":gibbs arithmetic is AbstractQAtlas's FreeEnergyLegendre" begin
    for (e, s, beta) in ((-1.3, 0.42, 0.5), (-0.75, 0.1, 2.0), (0.0, 0.0, 1.0))
        f_closed = e - s / beta   # the formula this edge used to carry inline
        f_derived = QAtlas.solve(QAtlas.FreeEnergyLegendre(), Val(:F); U=e, S=s, β=beta)
        @test f_derived ≈ f_closed atol = 1e-14
    end
    # β-or-T keyword convention must agree (T = 2 ⇔ β = 0.5).
    @test QAtlas.solve(QAtlas.FreeEnergyLegendre(), Val(:F); U=-1.3, S=0.42, T=2.0) ≈
        QAtlas.solve(QAtlas.FreeEnergyLegendre(), Val(:F); U=-1.3, S=0.42, β=0.5)
    # Exact inputs stay exact — the relation is algebraic, not floating-point.
    @test QAtlas.solve(QAtlas.FreeEnergyLegendre(), Val(:F); U=-1//2, S=1//4, β=1//2) ==
        -1//2 - (1//4) / (1//2)
end
