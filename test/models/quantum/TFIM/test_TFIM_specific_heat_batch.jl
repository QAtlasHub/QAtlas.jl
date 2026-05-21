# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_specific_heat_batch.jl
#
# TFIM specific heat per site at trivial limits:
#   * β → ∞ (T → 0): c = 0 exactly (gapped, exponentially small)
#   * β → 0 (T → ∞): c → 0 as ~ β² (high-T tail of finite-spectrum system)
# Pure verify(); branches off main. Refs #381.
#
# Tolerance notes:
#   * High-T tolerance 1e-4 is sized for c(β) ≈ β² · (J² + h²) / O(1).
#     At β=1e-3 and (J,h) up to (1,2): c ≈ (1e-3)² · 5 ≈ 5e-6, well inside 1e-4.
#   * Low-T β=1e2 (instead of 1e6) is sufficient: with min gap Δ ≈ 2|h−J| = 1,
#     β·Δ = 100 gives c ∝ (βΔ)²·exp(−βΔ) < 1e-30, far below 1e-9.
#
# PBC restriction:
#   PBC cards are restricted to h < J because test_TFIM_pbc_thermal.jl on main
#   documents an open parity-sector bug in the disordered phase (h > J). Running
#   trivial-limit cards in the broken regime would silently pass (c = 0 from
#   both sides) and defeat the purpose of verification.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — SpecificHeat trivial limits (#381 batch)" begin
    LOW_T_BETA = 1e2
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
            refs=[
                "TFIM at T → ∞: c → 0 as ~β² (high-T tail of bounded-spectrum quantum system)",
            ],
            fetch_kw=(; beta=HIGH_T_BETA),
        )
        for N in (8, 12)
            # Ordered-phase (h<J) finite-N OBC has a Z2-doublet splitting
            # delta ~ h^N (cat-state level pair); at N=8, h=0.5, delta ~ 0.004,
            # so beta*delta ~ 0.4 at LOW_T_BETA=100 and c is NOT exponentially
            # small. Restrict low-T finite-N c=0 claim to the disordered phase
            # (h>J) where Delta = 2(h-J) gives beta*Delta >> 1 truly.
            if h > J
                verify(
                    TFIM(; J=J, h=h),
                    SpecificHeat(),
                    OBC(N);
                    route=:limiting_case,
                    independent=0.0,
                    agree_within=1e-9,
                    refs=[
                        "TFIM OBC disordered phase (h>J), T → 0: c = 0 exactly (field-induced gap exponentially suppresses c)",
                    ],
                    fetch_kw=(; beta=LOW_T_BETA),
                )
            end
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
            # PBC finite-N low-T card dropped entirely:
            #   * h<J (ordered) finite-N PBC has Z2-doublet residue (same as OBC).
            #   * h>J (disordered) finite-N PBC has the documented parity-sector
            #     convention bug from test_TFIM_pbc_thermal.jl.
            # No PBC regime admits the c=0 trivial-limit claim at finite N here.
            # High-T (T -> infty) PBC card retained — finite N is fine in that
            # limit because the spectrum decouples from sector mixing.
            if h < J
                verify(
                    TFIM(; J=J, h=h),
                    SpecificHeat(),
                    PBC(N);
                    route=:limiting_case,
                    independent=0.0,
                    agree_within=1e-4,
                    refs=[
                        "TFIM PBC, T → ∞: c → 0 as β² (h < J regime to avoid parity-sector bug)",
                    ],
                    fetch_kw=(; beta=HIGH_T_BETA),
                )
            end
        end
    end
end
