# Shared runner for the generated-check suites (not a test_*.jl file — it is
# include()d by the per-slice files, so the universe completeness guard and
# the shard planner never schedule it directly).
#
# Contract: every emitted check passes; declared exclusions surface as
# visible :skip with their reason (no silent caps).

using QAtlas, Test
using QAtlas: run_generated_check

function run_generated_suite(checks; label::AbstractString)
    n_pass = n_skip = 0
    for c in checks
        out = run_generated_check(c)
        if out.status === :skip
            n_skip += 1
            @info "declared skip: $(c.id) — $(out.detail)"
            continue
        end
        out.status === :pass && (n_pass += 1)
        @testset "$(c.id)" begin
            if out.status !== :pass
                @info "generated check failed" c.id c.description out.lhs out.rhs out.abs_err out.rel_err out.detail
            end
            @test out.status === :pass
        end
    end
    println(
        label,
        ": ",
        n_pass,
        " pass / ",
        n_skip,
        " declared-skip / ",
        length(checks),
        " emitted",
    )
    return nothing
end
