# Generated derivation checks — the AbstractQAtlas relation network applied to
# what a hub fetches (src/core/derivation.jl).  Today's plane is `:scaling`: an
# exponent table held out one exponent at a time and solved for by every
# relation that reaches it from the rest.
#
# The anchors below are the reason this file is not just `run_generated_suite`.
# A generator whose routes all vanished would emit nothing and report a clean
# zero-failure run, so the surface is pinned: the relations that must appear,
# the hubs that must be reached, and the skips that must stay visible.

include("util_run_checks.jl")
using QAtlas: generated_checks, EXPONENT_SWEEPS, run_generated_check

@testset "generated derivation checks — :scaling" begin
    checks = generated_checks(; kinds=(:derivation,))
    @test !isempty(checks)
    ids = [c.id for c in checks]

    # Every SWEPT hub produces at least one real route, not merely an id under its
    # own name.  The bare prefix is satisfied by an excluded check too, so a hub
    # whose fetch started throwing would still match it — the id grammar
    # `<hub>/<point>/<target>/<relation>` is what separates a route from a refusal,
    # and a route is what this plane exists to emit.
    for s in EXPONENT_SWEEPS
        s isa QAtlas.SweptExponents || continue
        prefix = "derivation/scaling/$(QAtlas._kgshort(s.model))/"
        routes = filter(i -> startswith(i, prefix) && count(==('/'), i) >= 5, ids)
        @test !isempty(routes)
    end

    # All four `:scaling` routes fire somewhere.  A relation that dropped out
    # upstream (renamed, re-domained, made non-affine) would otherwise cost
    # coverage without costing a test.
    for r in ("Rushbrooke", "Widom", "Fisher", "Josephson")
        @test any(endswith("/$(r)"), ids)
    end

    # Both kinds of hub: the exact rational tables and the quoted-error ones.  Each
    # prefix runs to a sweep point, so only a real route matches it.
    for prefix in (
        "derivation/scaling/Universality{:Ising}/Infinite/d=2/",
        "derivation/scaling/Universality{:Ising}/Infinite/d=3/",
        "derivation/scaling/Universality{:Percolation}/Infinite/d=2/",
        "derivation/scaling/Universality{:Percolation}/Infinite/d=3/",
        "derivation/scaling/Universality{:Potts3}/Infinite/d=2/",
        "derivation/scaling/Universality{:MeanField}/Infinite/",
    )
        @test any(startswith(prefix), ids)
    end

    # `run_generated_suite` runs no assertion on a `:skip`, so a check that became
    # a skip by accident leaves the pass/fail tally untouched.  The skipped set is
    # therefore pinned whole rather than sampled: every skip is declared here, and
    # a new one fails this test instead of quietly removing a route.
    skipped = sort([c.id for c in checks if QAtlas.run_generated_check(c).status === :skip])
    @test skipped == sort([
        "derivation/scaling/Ising2D/Infinite",
        "derivation/scaling/KPZ1D/Infinite",
        "derivation/scaling/MeanField/Infinite",
        "derivation/scaling/Universality{:IsingSDRG}/Infinite",
        "derivation/scaling/Universality{:Ising}/Infinite/d=2/δ/Widom",
        "derivation/scaling/Universality{:Percolation}/Infinite/d=2/α/Rushbrooke",
        "derivation/scaling/Universality{:Percolation}/Infinite/d=2/γ/Fisher",
        "derivation/scaling/Universality{:Percolation}/Infinite/d=2/δ/Widom",
        "derivation/scaling/Universality{:XY}/Infinite/d=2",
    ])

    run_generated_suite(checks; label="generated derivation checks")
end

# The circular routes are the ones a green run would otherwise be built on: a
# value the table obtained FROM a relation cannot be checked BY it.  Assert they
# are skipped rather than passing, and that the reason travels with the skip.
@testset "a value derived from a relation is not checked by it" begin
    checks = generated_checks(; kinds=(:derivation,))
    by_id = Dict(c.id => c for c in checks)

    for id in (
        "derivation/scaling/Universality{:Ising}/Infinite/d=2/δ/Widom",
        "derivation/scaling/Universality{:Percolation}/Infinite/d=2/γ/Fisher",
        "derivation/scaling/Universality{:Percolation}/Infinite/d=2/δ/Widom",
        "derivation/scaling/Universality{:Percolation}/Infinite/d=2/α/Rushbrooke",
    )
        @test haskey(by_id, id)
        out = run_generated_check(by_id[id])
        @test out.status === :skip
        @test occursin("was obtained from", out.detail)
    end

    # ...and the OTHER routes to the same target still judge it, so declaring a
    # value derived costs one route, not the row.  Percolation d=2 reaches γ by
    # Rushbrooke and Widom besides the excluded Fisher.
    others = filter(
        c ->
            startswith(
                c.id, "derivation/scaling/Universality{:Percolation}/Infinite/d=2/γ/"
            ) && !endswith(c.id, "/Fisher"),
        checks,
    )
    @test !isempty(others)
    @test all(run_generated_check(c).status === :pass for c in others)
end

# A table whose letters mean something else must be kept OUT, and the refusal
# must be visible in the stream rather than an absence.
@testset "a table with different quantities behind the same letters is refused" begin
    checks = generated_checks(; kinds=(:derivation,))
    kpz = filter(c -> startswith(c.id, "derivation/scaling/KPZ1D/"), checks)
    @test length(kpz) == 1
    out = run_generated_check(only(kpz))
    @test out.status === :skip
    @test occursin("roughness exponent", out.detail)

    sdrg = filter(
        c -> startswith(c.id, "derivation/scaling/Universality{:IsingSDRG}/"), checks
    )
    @test length(sdrg) == 1
    @test run_generated_check(only(sdrg)).status === :skip
end

# The BKT point returns η alone.  Too few exponents to close is a REPORTED skip,
# not an empty result — an absent hub and a hub with nothing to check read the
# same in a pass count and must not in the stream.
@testset "a table too small to close is skipped, not absent" begin
    checks = generated_checks(; kinds=(:derivation,))
    bkt = filter(c -> c.id == "derivation/scaling/Universality{:XY}/Infinite/d=2", checks)
    @test length(bkt) == 1
    out = run_generated_check(only(bkt))
    @test out.status === :skip
    @test occursin("too few to close", out.detail)
end
