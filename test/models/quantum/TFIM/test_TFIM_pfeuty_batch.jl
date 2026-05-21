# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_pfeuty_batch.jl
#
# Pfeuty (1970) closed-form verification cards covering uncorroborated
# TFIM hubs. Pure verify() — no legacy here; self-contained, branches
# off main directly. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — Pfeuty closed-form cards (#381 batch)" begin
    # SpontaneousMagnetization/Infinite, ordered phase h < J:
    #   m_z(J,h) = (1 - (h/J)^2)^{1/8}  (Pfeuty 1970)
    # Points cover the classical-Ising limit (h=0 -> m_z=1), a mid-ordered
    # point, a J-scaled point, and a near-QCP point (h/J=0.95) that
    # stresses the 1/8 critical exponent where dm/dh diverges.
    for (J, h) in ((1.0, 0.0), (1.0, 0.5), (2.0, 1.0), (1.0, 0.95))
        verify(
            TFIM(; J=J, h=h),
            SpontaneousMagnetization(),
            Infinite();
            route=:second_closed_form,
            independent=(1 - (h / J)^2)^(1 / 8),
            agree_within=1e-12,
            refs=["Pfeuty 1970: m_z = (1 - (h/J)^2)^{1/8} for h < J (ordered phase)"],
        )
    end

    # SpontaneousMagnetization/Infinite, disordered / critical h ≥ J:
    #   m_z = 0  (Z₂ symmetry unbroken).
    # The QCP point (J=1, h=1) is intentionally included: at h=J both the
    # ordered-phase formula (1-(h/J)^2)^{1/8} and the disordered value 0
    # collapse to 0 exactly, so the card tests that src returns exactly 0
    # at the branch point (any non-zero residue would indicate a precision
    # or branch-selection bug in src).
    for (J, h) in ((1.0, 1.0), (1.0, 2.0), (0.5, 1.0))
        verify(
            TFIM(; J=J, h=h),
            SpontaneousMagnetization(),
            Infinite();
            route=:second_closed_form,
            independent=0.0,
            agree_within=1e-12,
            refs=["Pfeuty 1970: m_z = 0 for h ≥ J (disordered / critical)"],
        )
    end

    # MassGap/Infinite: Δ = 2|h − J| (Pfeuty 1970, thermodynamic-limit
    # quasi-particle gap from the BdG dispersion). The closed form ONLY
    # applies to the Infinite() dispatch — PBC(N) and OBC(N) drop back
    # to finite-N BdG / sector-comparison kernels whose values differ
    # substantially at small N (the ordered-phase doublet splitting is
    # exponentially small in N, and the critical-point gap is the CFT
    # 1/N correction). Card asserts the Infinite closed form.
    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0))
        verify(
            TFIM(; J=J, h=h),
            MassGap(),
            Infinite();
            route=:second_closed_form,
            independent=2 * abs(h - J),
            agree_within=1e-12,
            refs=[
                "Pfeuty 1970: Δ = 2|h − J| (thermodynamic-limit BdG gap; PBC/OBC kernels at finite N differ — see TFIM.jl docstring)",
            ],
        )
    end
end
