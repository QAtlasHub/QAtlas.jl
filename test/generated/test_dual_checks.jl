# Generated duality cross-checks (#699) — independent-implementation
# comparisons materialized from the @dual edges: fetch on the source, map
# parameters, fetch on the target, compare (through the edge's per-quantity
# value_map).  Both endpoint rows are independent (non-delegating) by
# construction — the structural antidote to circular verification.

using QAtlas, Test
using QAtlas: generated_checks, run_generated_check

@testset "generated duality cross-checks" begin
    checks = generated_checks(; kinds=(:dual,))
    @test !isempty(checks)

    ids = [c.id for c in checks]
    # the shipped catalog: TFIM Kramers–Wannier self-duality and the
    # Jordan–Wigner TFIM ↔ Kitaev wire map both materialize
    @test any(startswith("dual/tfim_kramers_wannier/"), ids)
    @test any(startswith("dual/tfim_kitaev_jordan_wigner/"), ids)

    n_pass = 0
    for c in checks
        out = run_generated_check(c)
        out.status === :pass && (n_pass += 1)
        @testset "$(c.id)" begin
            if out.status !== :pass
                @info "generated duality check failed" c.id c.description out.lhs out.rhs out.abs_err out.rel_err out.detail
            end
            @test out.status === :pass
        end
    end
    println("generated duality checks: $(n_pass) pass / $(length(checks)) emitted")
end
