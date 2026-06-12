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
