# Generated duality cross-checks (#699) — independent-implementation
# comparisons materialized from the @dual edges: fetch on the source, map
# parameters, fetch on the target, compare (through the edge's per-quantity
# value_map).  Both endpoint rows are independent (non-delegating) by
# construction — the structural antidote to circular verification.

include("util_run_checks.jl")
using QAtlas: generated_checks

@testset "generated duality cross-checks" begin
    checks = generated_checks(; kinds=(:dual,))
    @test !isempty(checks)

    ids = [c.id for c in checks]
    # the shipped catalog: TFIM Kramers–Wannier self-duality and the
    # Jordan–Wigner TFIM ↔ Kitaev wire map both materialize
    @test any(startswith("dual/tfim_kramers_wannier/"), ids)
    @test any(startswith("dual/tfim_kitaev_jordan_wigner/"), ids)

    run_generated_suite(checks; label="generated duality checks")
end
