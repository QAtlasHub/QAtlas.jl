# Generated identity checks not claimed by a dedicated slice file: a NEW
# @identity edge runs HERE automatically (no test edit needed) until it earns
# its own file.  The claimed-prefix list below must mirror
# test_identity_gibbs.jl and test_identity_isotropy.jl; an empty remainder is
# legal today.

include("util_run_checks.jl")
using QAtlas: generated_checks

const _CLAIMED_PREFIXES = ("identity/gibbs/", "identity/su2_")

@testset "generated identity checks — remainder" begin
    all_checks = generated_checks(; kinds=(:identity,))
    checks = filter(c -> !any(p -> startswith(c.id, p), _CLAIMED_PREFIXES), all_checks)
    # Each declared prefix must still match ≥ 1 check — a stale prefix (a slice
    # whose edge was renamed/removed) would otherwise silently route its
    # would-be checks here, hiding the drift.  (Tautological partition asserts
    # nothing; this catches the real failure mode.)
    for p in _CLAIMED_PREFIXES
        @test any(c -> startswith(c.id, p), all_checks)
    end
    run_generated_suite(checks; label="generated identity remainder")
end
