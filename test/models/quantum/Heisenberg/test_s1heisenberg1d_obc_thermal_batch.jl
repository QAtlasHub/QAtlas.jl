# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_s1heisenberg1d_obc_thermal_batch.jl
#
# Trivial-temperature-limit verification cards for the OBC S=1 Heisenberg
# chain. log 3 per site is the spin-1 paramagnet entropy density.
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "S1Heisenberg1D — OBC thermal trivial limits (#381 batch)" begin
    LOW_T_BETA = 1e6
    HIGH_T_BETA = 1e-3

    for J in (0.5, 1.0, 2.0)
        for N in (3, 4, 5)
            verify(
                S1Heisenberg1D(; J=J),
                ThermalEntropy(),
                OBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-9,
                refs=[
                    "S1Heisenberg1D OBC T → 0: unique GS at finite N ⇒ s = 0 exactly (Haldane gap controls the rate of approach in the thermodynamic limit)",
                ],
                fetch_kw=(; beta=LOW_T_BETA),
            )
            verify(
                S1Heisenberg1D(; J=J),
                ThermalEntropy(),
                OBC(N);
                route=:limiting_case,
                independent=log(3),
                agree_within=1e-5,
                refs=["S1Heisenberg1D OBC T → ∞: spin-1 paramagnet ⇒ s = log 3 per spin"],
                fetch_kw=(; beta=HIGH_T_BETA),
            )
            verify(
                S1Heisenberg1D(; J=J),
                SpecificHeat(),
                OBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-9,
                refs=["S1Heisenberg1D OBC T → 0: gap suppression ⇒ c = 0 exactly"],
                fetch_kw=(; beta=LOW_T_BETA),
            )
            verify(
                S1Heisenberg1D(; J=J),
                SpecificHeat(),
                OBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-4,
                refs=["S1Heisenberg1D OBC T → ∞: c → 0 as ~β² high-T tail"],
                fetch_kw=(; beta=HIGH_T_BETA),
            )
            verify(
                S1Heisenberg1D(; J=J),
                FreeEnergy(),
                OBC(N);
                route=:limiting_case,
                independent=-log(3) / HIGH_T_BETA,
                agree_within=1e-2,
                refs=[
                    "S1Heisenberg1D OBC T → ∞: spin-1 paramagnet f/N = -T log 3 = -log(3)/β"
                ],
                fetch_kw=(; beta=HIGH_T_BETA),
            )
        end
    end
end
