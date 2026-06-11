# Generated convergence-sequence checks (#701) — the executable form of the
# @limits_to edges: the source model, driven along the declared approach
# sequence, must converge to the (independently implemented) target value —
# shrinking error along the sequence, terminal error below the edge's
# final_atol.  The hand-written Δ → 1⁺ pattern of #691, as a generator.

using QAtlas, Test
using QAtlas: generated_checks, run_generated_check

@testset "generated limit convergence checks" begin
    checks = generated_checks(; kinds=(:limit,))
    @test !isempty(checks)

    ids = [c.id for c in checks]
    @test any(startswith("limit/xxz_isotropic_limit/"), ids)

    n_pass = 0
    for c in checks
        out = run_generated_check(c)
        out.status === :pass && (n_pass += 1)
        @testset "$(c.id)" begin
            if out.status !== :pass
                @info "generated limit check failed" c.id c.description out.lhs out.rhs out.abs_err out.detail
            end
            @test out.status === :pass
        end
    end
    println("generated limit checks: $(n_pass) pass / $(length(checks)) emitted")
end
