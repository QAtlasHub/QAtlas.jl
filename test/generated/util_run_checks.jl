# Shared runner for the generated-check suites (not a test_*.jl file — it is
# include()d by the per-slice files, so the universe completeness guard and
# the shard planner never schedule it directly).
#
# Contract: every emitted check passes; declared exclusions surface as
# visible :skip with their reason (no silent caps); a runner that throws is an
# :error (a config/dispatch bug, reported distinctly from a numerical :fail).

using QAtlas, Test
using QAtlas: run_generated_check

function run_generated_suite(checks; label::AbstractString)
    n_pass = n_skip = n_fail = n_error = 0
    for c in checks
        out = run_generated_check(c)
        if out.status === :skip
            n_skip += 1
            @info "declared skip: $(c.id) — $(out.detail)"
            continue
        end
        out.status === :pass && (n_pass += 1)
        out.status === :fail && (n_fail += 1)
        out.status === :error && (n_error += 1)
        @testset "$(c.id)" begin
            if out.status !== :pass
                # :fail = physics disagrees, :error = the runner threw (bug)
                @info "generated check $(out.status)" c.id c.description out.lhs out.rhs out.abs_err out.rel_err out.detail
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
        n_fail,
        " fail / ",
        n_error,
        " error / ",
        length(checks),
        " emitted",
    )
    return nothing
end
