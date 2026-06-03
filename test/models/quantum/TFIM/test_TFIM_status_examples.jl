# test/models/quantum/TFIM/test_TFIM_status_examples.jl
#
# Worked example for the registry `status` axis (v0.24):
#   * a one-sided :bound — the Lieb-Robinson velocity cone, saturated by
#     the maximum group velocity of the bare TFIM dispersion.
# Drives verify_bound against an INDEPENDENT witness (the bare dispersion),
# never against the same closed form.  (The :approx worked example is
# deferred to the definition-list redesign.)

using QAtlas, Test
using QAtlas: TFIM, LiebRobinsonBound, Infinite

@testset "TFIM LiebRobinsonBound — group velocity saturates v_LR (:bound)" begin
    # Dispersion Λ(k) = 2√(J² + h² − 2 J h cos k); group velocity is
    # v(k) = dΛ/dk = 2 J h sin k / √(J² + h² − 2 J h cos k).  Its maximum
    # over k is 2 min(|J|,|h|) (attained at cos k = min/max), exactly the
    # fetched Lieb-Robinson velocity.  So an independently computed max
    # group velocity stays ≤ v_LR and saturates it.
    for (J, h) in ((1.0, 0.5), (0.7, 1.3), (1.0, 1.0))
        m = TFIM(; J=J, h=h)
        vg(k) = abs(2 * J * h * sin(k)) / sqrt(J^2 + h^2 - 2 * J * h * cos(k) + 1e-300)
        ks = range(0, π; length=20001)
        max_vg = maximum(vg, ks)

        s = verify_bound(
            m,
            LiebRobinsonBound(),
            Infinite();
            route=:dispersion_velocity,
            measured=[max_vg],
            relation=:leq,
            saturating=true,
            slack=1e-3,
            refs=["LiebRobinson1972", "HastingsKoma2006"],
            at=["J=$J", "h=$h"],
        )
        @test s ≈ 2 * min(abs(J), abs(h)) atol = 1e-12   # subject = fetched v_LR
        @test max_vg <= s + 1e-9                          # cone genuinely bounds it
    end
end
