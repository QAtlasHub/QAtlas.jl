# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_xxcorrelation_inf_batch.jl
#
# Two trivial limits of the TFIM XXCorrelation in the thermodynamic limit:
#   * i = j: ⟨σ^x_i σ^x_i⟩ = ⟨I⟩ = 1 exactly (Pauli σ_x² = I), for ANY
#            (J, h, β).
#   * i ≠ j at T → ∞: the maximally mixed state ρ = I/2^N gives
#            ⟨σ^x_i σ^x_j⟩ = 0 for any i ≠ j (uncorrelated paramagnet).
# Pure verify(); branches off main. Refs #381.
#
# TODO(#438): this card uses the σ (Pauli) convention for XXCorrelation. Once
# PR #438 merges and the hub switches to the S (spin) convention, divide all
# `independent` targets by 4 (⟨S^α_i S^α_j⟩ = ⟨σ^α_i σ^α_j⟩ / 4 for spin-1/2):
# the i=j target becomes 1/4 and the high-T i≠j target stays 0.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — XXCorrelation/Infinite trivial limits (#381 batch)" begin
    # i = j: ⟨σ^x²⟩ = 1 for all (J, h, β)
    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (0.5, 2.0))
        for β in (1e-3, 10.0, 1e6)
            verify(
                TFIM(; J=J, h=h),
                XXCorrelation(mode=:static),
                Infinite();
                route=:second_closed_form,
                independent=1.0,
                agree_within=1e-10,
                refs=["⟨σ^x_i σ^x_i⟩ = ⟨I⟩ = 1 exactly (σ_x² = I), independent of state"],
                fetch_kw=(; beta=β, i=40, j=40),
            )
        end
    end

    # i ≠ j at T → ∞: uncorrelated paramagnet ⇒ ⟨σ^x_i σ^x_j⟩ → 0
    HIGH_T_BETA = 1e-3
    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0))
        for (i, j) in ((40, 41), (40, 50), (40, 60))
            verify(
                TFIM(; J=J, h=h),
                XXCorrelation(mode=:static),
                Infinite();
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-3,
                refs=["TFIM at T → ∞: ρ = I/2^N ⇒ ⟨σ^x_i σ^x_j⟩ = 0 for i ≠ j (uncorrelated paramagnet)"],
                fetch_kw=(; beta=HIGH_T_BETA, i=i, j=j),
            )
        end
    end
end
