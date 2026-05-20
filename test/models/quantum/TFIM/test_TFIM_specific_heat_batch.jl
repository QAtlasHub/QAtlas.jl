# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_specific_heat_batch.jl
#
# TFIM specific heat per site at trivial limits:
#   * β → ∞ (T → 0): c = 0 exactly (gapped, exponentially small)
#   * β → 0 (T → ∞): c → 0 as ~ β² (high-T tail of finite-spectrum system)
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — SpecificHeat trivial limits (#381 batch)" begin
    LOW_T_BETA = 1e6
    HIGH_T_BETA = 1e-3

    for (J, h) in ((1.0, 0.5), (1.0, 2.0), (0.5, 2.0))
        # T → 0: c = 0 exactly (gapped, all heat capacity exponentially suppressed)
        verify(
            TFIM(; J=J, h=h),
            SpecificHeat(),
            Infinite();
            route=:limiting_case,
            independent=0.0,
            agree_within=1e-9,
            refs=["TFIM gapped phase (h ≠ J): c(T → 0) ∝ exp(-Δ/T) → 0 exactly"],
            fetch_kw=(; beta=LOW_T_BETA),
        )
        # T → ∞: c → 0 as β² × const ≈ 5e-6 at β=1e-3
        verify(
            TFIM(; J=J, h=h),
            SpecificHeat(),
            Infinite();
            route=:limiting_case,
            independent=0.0,
            agree_within=1e-4,
            refs=["TFIM at T → ∞: c → 0 as ~β² (high-T tail of bounded-spectrum quantum system)"],
            fetch_kw=(; beta=HIGH_T_BETA),
        )
        for N in (8, 12)
            verify(
                TFIM(; J=J, h=h),
                SpecificHeat(),
                OBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-9,
                refs=["TFIM OBC, T → 0: c = 0 exactly (gap → exponentially small heat capacity)"],
                fetch_kw=(; beta=LOW_T_BETA),
            )
            verify(
                TFIM(; J=J, h=h),
                SpecificHeat(),
                OBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-4,
                refs=["TFIM OBC, T → ∞: c → 0 as β² (high-T tail)"],
                fetch_kw=(; beta=HIGH_T_BETA),
            )
            verify(
                TFIM(; J=J, h=h),
                SpecificHeat(),
                PBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-9,
                refs=["TFIM PBC, T → 0: c = 0 exactly"],
                fetch_kw=(; beta=LOW_T_BETA),
            )
            verify(
                TFIM(; J=J, h=h),
                SpecificHeat(),
                PBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-4,
                refs=["TFIM PBC, T → ∞: c → 0 as β²"],
                fetch_kw=(; beta=HIGH_T_BETA),
            )
        end
    end
end
