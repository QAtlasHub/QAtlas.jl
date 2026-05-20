# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_thermal_entropy_batch.jl
#
# TFIM thermal entropy per site at trivial limits:
#   * β → ∞ (T → 0): s = 0 exactly (unique gapped GS, no residual entropy)
#   * β → 0 (T → ∞): s → log 2 (paramagnet, 2 states per site)
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — ThermalEntropy trivial limits (#381 batch)" begin
    LOW_T_BETA = 1e6   # β → ∞
    HIGH_T_BETA = 1e-3 # β → 0

    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (0.5, 1.0))
        # /Infinite, OBC, PBC at T → 0: s = 0 exactly (unique GS, gapped at h ≠ J)
        for N in (8, 12, 16)
            verify(
                TFIM(; J=J, h=h),
                ThermalEntropy(),
                Infinite();
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-9,
                refs=["TFIM is gapped for h ≠ J ⇒ unique GS ⇒ thermal entropy density vanishes as T → 0"],
                fetch_kw=(; beta=LOW_T_BETA),
            )
            verify(
                TFIM(; J=J, h=h),
                ThermalEntropy(),
                OBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-9,
                refs=["TFIM OBC unique GS at T → 0 ⇒ s = 0 exactly"],
                fetch_kw=(; beta=LOW_T_BETA),
            )
            verify(
                TFIM(; J=J, h=h),
                ThermalEntropy(),
                PBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-9,
                refs=["TFIM PBC unique GS at T → 0 ⇒ s = 0 exactly"],
                fetch_kw=(; beta=LOW_T_BETA),
            )
        end
    end

    # High-T limit: s → log 2 per site (paramagnet, 2 states/spin).
    for (J, h) in ((1.0, 1.0),)
        verify(
            TFIM(; J=J, h=h),
            ThermalEntropy(),
            Infinite();
            route=:limiting_case,
            independent=log(2),
            agree_within=1e-5,
            refs=["TFIM at T → ∞: s → log 2 per site (paramagnet)"],
            fetch_kw=(; beta=HIGH_T_BETA),
        )
        for N in (8, 12, 16)
            verify(
                TFIM(; J=J, h=h),
                ThermalEntropy(),
                OBC(N);
                route=:limiting_case,
                independent=log(2),
                agree_within=1e-5,
                refs=["TFIM OBC at T → ∞: s → log 2 per site"],
                fetch_kw=(; beta=HIGH_T_BETA),
            )
            verify(
                TFIM(; J=J, h=h),
                ThermalEntropy(),
                PBC(N);
                route=:limiting_case,
                independent=log(2),
                agree_within=1e-5,
                refs=["TFIM PBC at T → ∞: s → log 2 per site"],
                fetch_kw=(; beta=HIGH_T_BETA),
            )
        end
    end
end
