# Generated symmetry corroboration checks (#700) — the :symmetry kind: every
# @symmetry profile declaring a `gapped` fact is cross-checked against the
# model's registered MassGap implementation at Infinite, closing the
# two-sources-of-truth seam between the profile store and REGISTRY.

include("util_run_checks.jl")
using QAtlas: generated_checks

@testset "generated symmetry corroboration checks" begin
    checks = generated_checks(; kinds=(:symmetry,))
    @test !isempty(checks)

    ids = [c.id for c in checks]
    # both directions of the claim are exercised: a declared-gapless profile
    # (Heisenberg1D) and a declared-gapped one (S1Heisenberg1D Haldane gap)
    @test any(startswith("symmetry/gapped/Heisenberg1D/"), ids)
    @test any(startswith("symmetry/gapped/S1Heisenberg1D/"), ids)

    run_generated_suite(checks; label="generated symmetry checks")
end
