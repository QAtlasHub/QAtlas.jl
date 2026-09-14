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

    # Every declared hub reaches the check stream under its own name, swept or
    # refused.  A sweep that silently stopped generating is the failure mode a
    # zero-failure run cannot distinguish from success.
    for s in EXPONENT_SWEEPS
        prefix = "derivation/scaling/$(QAtlas._kgshort(s.model))/"
        @test any(startswith(prefix), ids)
    end

    # All four `:scaling` routes fire somewhere.  A relation that dropped out
    # upstream (renamed, re-domained, made non-affine) would otherwise cost
    # coverage without costing a test.
    for r in ("Rushbrooke", "Widom", "Fisher", "Josephson")
        @test any(endswith("/$(r)"), ids)
    end

    # Both kinds of hub: the exact rational tables and the quoted-error ones.
    @test any(startswith("derivation/scaling/Universality{:Ising}/Infinite/d=2/"), ids)
    @test any(startswith("derivation/scaling/Universality{:Ising}/Infinite/d=3/"), ids)
    @test any(startswith("derivation/scaling/Universality{:Percolation}/"), ids)
    @test any(startswith("derivation/scaling/Universality{:Potts3}/"), ids)
    @test any(startswith("derivation/scaling/Universality{:MeanField}/"), ids)

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
