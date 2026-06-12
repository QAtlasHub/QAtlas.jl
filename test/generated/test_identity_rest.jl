# Generated identity checks not claimed by a dedicated slice file: a NEW
# @identity edge runs HERE automatically (no test edit needed) until it earns
# its own file.  The claimed-prefix list below must mirror
# test_identity_gibbs.jl and test_identity_isotropy.jl; an empty remainder is
# legal today.

include("util_run_checks.jl")
using QAtlas: generated_checks

const _CLAIMED_PREFIXES = ("identity/gibbs/", "identity/su2_")

@testset "generated identity checks — remainder" begin
    all_ids = generated_checks(; kinds=(:identity,))
    checks = filter(c -> !any(p -> startswith(c.id, p), _CLAIMED_PREFIXES), all_ids)
    # every identity check is claimed by exactly one slice file
    @test length(checks) +
          count(c -> any(p -> startswith(c.id, p), _CLAIMED_PREFIXES), all_ids) ==
        length(all_ids)
    run_generated_suite(checks; label="generated identity remainder")
end
