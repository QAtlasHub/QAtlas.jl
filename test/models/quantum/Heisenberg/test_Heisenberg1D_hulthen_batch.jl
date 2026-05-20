# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_Heisenberg1D_hulthen_batch.jl
#
# Hulthén (1938) Bethe-ansatz closed-form verification cards for the
# spin-½ Heisenberg chain. Pure verify(); self-contained, branches off
# main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Heisenberg1D — Hulthén closed-form cards (#381 batch)" begin
    # GroundStateEnergyDensity/Infinite: e₀ = J(1/4 − log 2) (Hulthén 1938,
    # Bethe-ansatz exact). Linear J-scaling.
    for J in (0.5, 1.0, 2.0)
        verify(
            Heisenberg1D(; J=J),
            GroundStateEnergyDensity(),
            Infinite();
            route=:second_closed_form,
            independent=J * (1 / 4 - log(2)),
            agree_within=1e-14,
            refs=["Hulthén 1938: e₀ = J(1/4 − log 2) Bethe-ansatz exact"],
        )
    end
end
