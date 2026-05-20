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
    for (J, h) in ((1.0, 0.0), (1.0, 0.5), (2.0, 1.0))
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

    # MassGap/PBC: Δ = 2|h − J|, BC-independent (Pfeuty 1970).
    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0))
        verify(
            TFIM(; J=J, h=h),
            MassGap(),
            PBC(8);
            route=:second_closed_form,
            independent=2 * abs(h - J),
            agree_within=1e-12,
            refs=["Pfeuty 1970: Δ = 2|h − J| (BC-independent)"],
        )
    end
end
