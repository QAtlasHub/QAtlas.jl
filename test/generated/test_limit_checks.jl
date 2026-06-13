# Generated convergence-sequence checks (#701) — the executable form of the
# @limits_to edges: the source model, driven along the declared approach
# sequence, must converge to the (independently implemented) target value —
# shrinking error along the sequence (within the edge's declared mono_slack),
# terminal error below the spec's final_atol.  The hand-written Δ-sequence
# pattern of #691, as a generator.

include("util_run_checks.jl")
using QAtlas: generated_checks

@testset "generated limit convergence checks" begin
    checks = generated_checks(; kinds=(:limit,))
    @test !isempty(checks)

    ids = [c.id for c in checks]
    @test any(startswith("limit/xxz_isotropic_limit/"), ids)

    run_generated_suite(checks; label="generated limit checks")
end
