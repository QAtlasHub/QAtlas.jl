# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/XXZ/test_xxz1d_obc_thermal_batch.jl
#
# Trivial-temperature-limit verification cards for the OBC XXZ1D chain
# (spin-1/2). Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "XXZ1D — OBC thermal trivial limits (#381 batch)" begin
    LOW_T_BETA = 1e6
    HIGH_T_BETA = 1e-3

    for (J, Δ) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0))
        for N in (4, 6, 8)
            verify(
                XXZ1D(),
                ThermalEntropy(),
                OBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-9,
                refs=["XXZ1D OBC T → 0: unique-singlet GS ⇒ s = 0 exactly"],
                fetch_kw=(; J=J, Δ=Δ, beta=LOW_T_BETA),
            )
            verify(
                XXZ1D(),
                ThermalEntropy(),
                OBC(N);
                route=:limiting_case,
                independent=log(2),
                agree_within=1e-5,
                refs=["XXZ1D OBC T → ∞: spin-1/2 paramagnet ⇒ s = log 2 per spin"],
                fetch_kw=(; J=J, Δ=Δ, beta=HIGH_T_BETA),
            )
            verify(
                XXZ1D(),
                SpecificHeat(),
                OBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-9,
                refs=["XXZ1D OBC T → 0: c = 0 exactly"],
                fetch_kw=(; J=J, Δ=Δ, beta=LOW_T_BETA),
            )
            verify(
                XXZ1D(),
                SpecificHeat(),
                OBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-4,
                refs=["XXZ1D OBC T → ∞: c → 0 as ~β² high-T tail"],
                fetch_kw=(; J=J, Δ=Δ, beta=HIGH_T_BETA),
            )
            verify(
                XXZ1D(),
                FreeEnergy(),
                OBC(N);
                route=:limiting_case,
                independent=-log(2) / HIGH_T_BETA,
                agree_within=1e-2,
                refs=["XXZ1D OBC T → ∞: free paramagnet f/N = -T log 2 = -log(2)/β"],
                fetch_kw=(; J=J, Δ=Δ, beta=HIGH_T_BETA),
            )
        end
    end
end
