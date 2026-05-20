# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_correlation_length_batch.jl
#
# Pfeuty/Bogoliubov closed-form verification card for the TFIM
# correlation length in the gapped phases (h ≠ J). Pure verify();
# self-contained, branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — CorrelationLength Pfeuty closed form (#381 batch)" begin
    # Pfeuty/Bogoliubov: ξ = 1 / log(max(J,h) / min(J,h)) for h ≠ J
    # (gapped ordered/disordered phases; TFIM self-duality gives the
    # same ξ at (J,h) and (h,J)). h = J is the critical point ξ = ∞,
    # excluded.
    for (J, h) in ((1.0, 0.5), (1.0, 2.0), (2.0, 1.0))
        ξ_closed = 1 / log(max(J, h) / min(J, h))
        verify(
            TFIM(; J=J, h=h),
            CorrelationLength(),
            Infinite();
            route=:second_closed_form,
            independent=ξ_closed,
            agree_within=1e-12,
            refs=["Pfeuty 1970: ξ = 1/log(max(J,h)/min(J,h)) (TFIM gapped phases)"],
        )
    end
end
