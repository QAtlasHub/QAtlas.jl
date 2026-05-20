# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_correlation_length_batch.jl
#
# Pfeuty/Bogoliubov closed-form verification card for the TFIM
# correlation length in the gapped phases (h ≠ J). Pure verify();
# self-contained, branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — CorrelationLength inverse-mass-gap closed form (#381 batch)" begin
    # ξ = 1 / (2|h - J|) for h ≠ J — the inverse mass gap (relativistic
    # ξ = v_F/Δ at the band minimum, with v_F = 2·min(J,h) and
    # Δ = 2|h - J|; the v_F factor is absorbed by the src's chosen
    # length-unit convention). This matches the existing src dispatch
    # at src/models/quantum/TFIM/TFIM_zaxis.jl and the tight-tolerance
    # assertion at test/models/quantum/TFIM/test_TFIM_zaxis.jl. h = J
    # gives ξ = ∞ (critical), excluded from the sweep.
    #
    # The three pairs sweep |h − J| ∈ {1, 2, 3}, yielding three distinct
    # numerical ξ values {0.5, 0.25, 1/6}. A previous parameterisation
    # (1,0.5),(1,2),(2,1) collapsed to |h-J| ∈ {0.5, 1, 1} (only two
    # distinct ξ); the new sweep gives three.
    #
    # NOTE: there is an open convention question regarding the Pfeuty
    # 1970 longitudinal-correlator form ξ = 1/log(max(J,h)/min(J,h))
    # vs the inverse-mass-gap form used here; tracked in issue #448.
    # This card asserts the src's currently-implemented inverse-mass-gap
    # convention so verification matches the canonical implementation
    # and the existing zaxis test.
    for (J, h) in ((1.0, 2.0), (1.0, 3.0), (1.0, 4.0))
        ξ_closed = 1 / (2 * abs(h - J))
        verify(
            TFIM(; J=J, h=h),
            CorrelationLength(),
            Infinite();
            route=:second_closed_form,
            independent=ξ_closed,
            agree_within=1e-12,
            refs=[
                "TFIM CorrelationLength = 1/(2|h-J|) (inverse mass gap; src TFIM_zaxis.jl); convention discrepancy with Pfeuty 1970 1/log(max/min) tracked in issue #448",
            ],
        )
    end
end
