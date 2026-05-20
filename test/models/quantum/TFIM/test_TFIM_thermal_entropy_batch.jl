# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_thermal_entropy_batch.jl
#
# TFIM thermal entropy per site at trivial limits:
#   * β → ∞ (T → 0): s = 0 exactly (unique gapped GS, no residual entropy)
#   * β → 0 (T → ∞): s → log 2 (paramagnet, 2 states per site)
# Pure verify(); branches off main. Refs #381.
#
# Scope choices (from review):
#   - Critical point (J=h) is excluded from the T→0 sweep — at the QCP the
#     gap closes and s does not vanish at finite β (CFT residual ~ π·c·v/(3β)).
#   - PBC cards are restricted to h < J. test_TFIM_pbc_thermal.jl documents a
#     parity-sector convention bug in the disordered phase h > J that needs a
#     separate sector-projection rewrite before PBC cards can include h ≥ J.
#   - LOW_T_BETA = 1e3 is large enough (β·Δ ≫ 1 for the gapped phases) and
#     avoids the numerical pathologies of β = 1e6 with non-shifted BdG sums.
#   - agree_within = 1e-5 for OBC/PBC to accommodate Z₂-doublet finite-size
#     splitting Δ ~ J(h/J)^N: at N=16, h=0.5 the doublet contributes ~4e-6
#     per site to s, which would violate a tighter 1e-9 tolerance.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — ThermalEntropy trivial limits (#381 batch)" begin
    LOW_T_BETA = 1e3   # β → ∞ (β·Δ ≫ 1 suffices for gapped phases)
    HIGH_T_BETA = 1e-3 # β → 0

    # T → 0 sweep: gapped (J, h) only — (1.0, 1.0) is the critical point and
    # is excluded. PBC additionally restricted to h < J (see header).
    for (J, h) in ((1.0, 0.5), (1.0, 2.0), (0.5, 1.0))
        for N in (8, 12, 16)
            verify(
                TFIM(; J=J, h=h),
                ThermalEntropy(),
                Infinite();
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-9,
                refs=[
                    "TFIM is gapped for h ≠ J ⇒ unique GS ⇒ thermal entropy density vanishes as T → 0",
                ],
                fetch_kw=(; beta=LOW_T_BETA),
            )
            verify(
                TFIM(; J=J, h=h),
                ThermalEntropy(),
                OBC(; N=N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-5,
                refs=[
                    "TFIM OBC unique GS at T → 0 ⇒ s = 0 (tolerance 1e-5 absorbs Z₂-doublet finite-size splitting)",
                ],
                fetch_kw=(; beta=LOW_T_BETA),
            )
            # PBC only in the ordered phase h < J (parity-sector bug in disordered phase).
            if h < J
                verify(
                    TFIM(; J=J, h=h),
                    ThermalEntropy(),
                    PBC(; N=N);
                    route=:limiting_case,
                    independent=0.0,
                    agree_within=1e-5,
                    refs=[
                        "TFIM PBC ordered phase (h<J) unique GS at T → 0 ⇒ s = 0 (tolerance 1e-5 absorbs Z₂-doublet splitting)",
                    ],
                    fetch_kw=(; beta=LOW_T_BETA),
                )
            end
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
                OBC(; N=N);
                route=:limiting_case,
                independent=log(2),
                agree_within=1e-5,
                refs=["TFIM OBC at T → ∞: s → log 2 per site"],
                fetch_kw=(; beta=HIGH_T_BETA),
            )
            verify(
                TFIM(; J=J, h=h),
                ThermalEntropy(),
                PBC(; N=N);
                route=:limiting_case,
                independent=log(2),
                agree_within=1e-5,
                refs=["TFIM PBC at T → ∞: s → log 2 per site"],
                fetch_kw=(; beta=HIGH_T_BETA),
            )
        end
    end
end
