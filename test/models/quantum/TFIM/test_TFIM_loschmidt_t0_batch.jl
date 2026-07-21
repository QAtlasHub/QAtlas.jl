# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_loschmidt_t0_batch.jl
#
# Loschmidt-echo rate function r(t) = -(1/N) log |⟨ψ₀|e^{-iHt}|ψ₀⟩|² at
# t = 0: the time-evolution operator is identity, so ⟨ψ₀|ψ₀⟩ = 1 and
# r(0) = 0 identically, for ANY initial-final TFIM pair.
# Follow-up: a non-trivial h_i ≠ h_f quench at finite t (probing DQPT
# zeros of the Loschmidt amplitude across the quantum critical point
# h_c = J) would be the natural next regression card.
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — LoschmidtEcho rate function at t=0 = 0 (#381 batch)" begin
    # /Infinite
    for (J, h_i, h_f) in (
        (1.0, 0.5, 1.0),
        (1.0, 0.5, 2.0),
        (1.0, 1.5, 1.0),
        (1.0, 1.5, 2.0),
        (2.0, 0.5, 3.0),
        (0.5, 2.0, 1.0),
    )
        verify(
            TFIM(; J=J, h=h_f),
            LoschmidtRateFunction(),
            Infinite();
            route=:second_closed_form,
            independent=0.0,
            agree_within=1e-10,
            refs=[
                "Loschmidt rate function r(t) = -(1/N) log |⟨ψ₀|e^{-iHt}|ψ₀⟩|² ⇒ r(t=0) = 0 (identity evolution, normalized state)",
            ],
            fetch_kw=(; initial=TFIM(; J=J, h=h_i), t=0.0),
        )
        # /OBC
        for N in (8, 12)
            verify(
                TFIM(; J=J, h=h_f),
                LoschmidtRateFunction(),
                OBC(N);
                route=:second_closed_form,
                independent=0.0,
                agree_within=1e-10,
                refs=[
                    "Loschmidt rate function r(t=0) = 0 in any BC, for any quench (initial state normalized)",
                ],
                fetch_kw=(; initial=TFIM(; J=J, h=h_i), t=0.0),
            )
        end
    end
end
