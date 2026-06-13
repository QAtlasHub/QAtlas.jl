# Generated SU(2) isotropy checks (#698 × #700 × #690) — the family-mode
# slice of the @identity edges (id prefix "identity/su2_"): pairwise
# component equalities over the taxonomy supertypes, gated on the @symmetry
# profiles.  Split from the Gibbs slice for shard-planner granularity; the
# prefix list must stay in sync with test_identity_rest.jl.

include("util_run_checks.jl")
using QAtlas: generated_checks

@testset "generated identity checks — SU(2) isotropy" begin
    checks = filter(
        c -> startswith(c.id, "identity/su2_"), generated_checks(; kinds=(:identity,))
    )
    @test !isempty(checks)

    # the @symmetry gate: SU(2)-profiled models generate, U(1) families do not
    ids = [c.id for c in checks]
    @test any(startswith("identity/su2_susceptibility_isotropy/Heisenberg1D/"), ids)
    @test !any(startswith("identity/su2_susceptibility_isotropy/XXZ1D/"), ids)

    run_generated_suite(checks; label="generated isotropy checks")
end
