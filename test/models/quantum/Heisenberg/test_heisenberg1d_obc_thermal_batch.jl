# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_heisenberg1d_obc_thermal_batch.jl
#
# Trivial-temperature-limit verification cards for the OBC Heisenberg
# chain:
#   ThermalEntropy : T → 0 ⇒ s = 0 (singlet GS); T → ∞ ⇒ s = log 2
#                    (paramagnet, 2 states/spin)
#   SpecificHeat   : T → 0 ⇒ c = 0 (gap suppression); T → ∞ ⇒ c → 0 (β² tail)
#   FreeEnergy     : T → ∞ ⇒ f/N → -T log 2 = -log(2)/β (paramagnet)
#
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Heisenberg1D — OBC thermal trivial limits (#381 batch)" begin
    LOW_T_BETA = 1e6
    HIGH_T_BETA = 1e-3

    for J in (0.5, 1.0, 2.0)
        for N in (4, 6, 8)
            # ThermalEntropy
            verify(
                Heisenberg1D(),
                ThermalEntropy(),
                OBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-9,
                refs=["Heisenberg1D OBC T → 0: spin-singlet GS ⇒ s = 0 exactly"],
                fetch_kw=(; J=J, beta=LOW_T_BETA),
            )
            verify(
                Heisenberg1D(),
                ThermalEntropy(),
                OBC(N);
                route=:limiting_case,
                independent=log(2),
                agree_within=1e-5,
                refs=["Heisenberg1D OBC T → ∞: free paramagnet ⇒ s = log 2 per spin"],
                fetch_kw=(; J=J, beta=HIGH_T_BETA),
            )
            # SpecificHeat
            verify(
                Heisenberg1D(),
                SpecificHeat(),
                OBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-9,
                refs=["Heisenberg1D OBC T → 0: gap suppression ⇒ c = 0 exactly"],
                fetch_kw=(; J=J, beta=LOW_T_BETA),
            )
            verify(
                Heisenberg1D(),
                SpecificHeat(),
                OBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-4,
                refs=["Heisenberg1D OBC T → ∞: c → 0 as ~β² (high-T tail)"],
                fetch_kw=(; J=J, beta=HIGH_T_BETA),
            )
            # FreeEnergy at T → ∞: f/N → -T log 2 = -log(2)/β
            verify(
                Heisenberg1D(),
                FreeEnergy(),
                OBC(N);
                route=:limiting_case,
                independent=-log(2) / HIGH_T_BETA,
                agree_within=1e-2,
                refs=["Heisenberg1D OBC T → ∞: free paramagnet ⇒ f/N = -T log 2 = -log(2)/β"],
                fetch_kw=(; J=J, beta=HIGH_T_BETA),
            )
        end
    end
end
