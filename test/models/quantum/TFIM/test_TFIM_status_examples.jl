# test/models/quantum/TFIM/test_TFIM_status_examples.jl
#
# Worked examples for the registry `status` axis (v0.24):
#   * a one-sided :bound  — the Lieb-Robinson velocity cone, saturated by
#     the maximum group velocity of the bare TFIM dispersion;
#   * a domain-limited :approx — the high-temperature (small-β) expansion
#     of the per-site free energy, tracked against the exact BdG value.
# Both drive verify_bound / verify_approx against an INDEPENDENT witness
# (dispersion / exact free energy), never against the same closed form.

using QAtlas, Test
using QAtlas: TFIM, LiebRobinsonBound, HighTemperatureFreeEnergy, FreeEnergy, Infinite

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

@testset "TFIM HighTemperatureFreeEnergy — small-β expansion (:approx)" begin
    # f(β)/N ≈ −ln2/β − (β/2)(J²+h²) + O(β³).  In-domain (βJ ≪ 1) it tracks
    # the exact BdG free energy; the reference is the exact fetch and the
    # tolerance reflects the O(β³) truncation.
    m = TFIM(; J=1.0, h=0.5)
    β = 0.05
    f_exact = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)

    s = verify_approx(
        m,
        HighTemperatureFreeEnergy(),
        Infinite();
        route=:high_temperature,
        reference=f_exact,
        agree_within=1e-3,
        valid_domain="betaJ << 1",
        error_order="O(beta^3)",
        refs=["Pfeuty1970"],
        fetch_kw=(; beta=β),
        at=["beta=$β"],
    )
    @test isapprox(s, f_exact; atol=1e-3)     # in-domain agreement

    # Out of domain (β ~ 1) the approximation must visibly break — this
    # guards against silently widening the stated validity window.
    f_exact_cold = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1.0)
    f_approx_cold = QAtlas.fetch(m, HighTemperatureFreeEnergy(), Infinite(); beta=1.0)
    @test !isapprox(f_approx_cold, f_exact_cold; atol=1e-3)
end
