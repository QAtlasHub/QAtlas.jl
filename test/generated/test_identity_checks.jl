# Generated identity checks (#698) — the executable form of the @identity
# edges, materialized by `generated_checks` against the live registry.
#
# Every (edge × hub × sweep point) the generator emits must pass; declared
# exclusions surface as :skip with their reason (visible accounting — no
# silent caps).  The hand-written harness (test/util/thermodynamic_identities
# .jl, exercised in test/identities/) stays the reference implementation
# until full parity is reached (#698's migration condition) — the anchors
# below pin the generated surface to cover at least the harness's Gibbs hubs.

using QAtlas, Test
using QAtlas: generated_checks, run_generated_check

@testset "generated identity checks" begin
    checks = generated_checks(; kinds=(:identity,))
    @test !isempty(checks)

    # Parity anchors: the generated Gibbs surface covers the hand-written
    # harness's hubs (TFIM at OBC + Infinite), and reaches hubs that have NO
    # hand-written identity file — the "new models get identity coverage for
    # free" payoff.
    ids = [c.id for c in checks]
    @test any(startswith("identity/gibbs/TFIM/Infinite"), ids)
    @test any(startswith("identity/gibbs/TFIM/OBC"), ids)
    @test any(startswith("identity/gibbs/XXZ1D/"), ids)
    # symmetry-gated isotropy fires only for @symmetry internal=:SU2 models
    @test any(startswith("identity/su2_susceptibility_isotropy/Heisenberg1D/"), ids)
    @test !any(startswith("identity/su2_susceptibility_isotropy/XXZ1D/"), ids)

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
                @info "generated identity check failed" c.id c.description out.lhs out.rhs out.abs_err out.rel_err out.detail
            end
            @test out.status === :pass
        end
    end
    println(
        "generated identity checks: $(n_pass) pass / $(n_skip) declared-skip / ",
        "$(length(checks)) emitted",
    )
end
