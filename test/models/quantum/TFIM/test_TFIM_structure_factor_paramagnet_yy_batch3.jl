# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_structure_factor_paramagnet_yy_batch3.jl
#
# Part 3/4 of the YY-axis paramagnet structure-factor batch (one (J, h) pair
# per file so the dense OBC(12) thermal ED balances finely across CI shards).
# Same paramagnet physics: at T → ∞, ρ = I / 2^N ⇒ ⟨σ^yy_i σ^yy_j⟩ = δ_ij and
# S_yyyy(q) = 1 for all q, all BC.
#
# TODO(#438): once #438 lands, migrate to S = σ/2 units (independent = 1/4,
# agree_within 1e-9). Until then σ convention is internally consistent.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — YYStructureFactor at T→∞ paramagnet = 1 (#381 batch, part 3/4)" begin
    HIGH_T_BETA = 1e-3

    for (J, h) in ((1.0, 2.0),)
        q = YYStructureFactor()
        for (q_name, q_val) in (("q=0", 0.0), ("q=π/2", π/2), ("q=π", π))
            # /Infinite
            verify(
                TFIM(; J=J, h=h),
                q,
                Infinite();
                route=:limiting_case,
                independent=1.0,
                agree_within=1e-2,
                refs=[
                    "TFIM at T → ∞: ρ = I/2^N ⇒ ⟨σ^α_i σ^α_j⟩ = δ_{ij} ⇒ S_αα(q) = 1 for any q, any axis",
                ],
                fetch_kw=(; q=q_val, beta=HIGH_T_BETA),
            )
            # /OBC (uses bc.N for the system size)
            for N in (8, 12)
                verify(
                    TFIM(; J=J, h=h),
                    q,
                    OBC(N);
                    route=:limiting_case,
                    independent=1.0,
                    agree_within=1e-2,
                    refs=[
                        "TFIM OBC at T → ∞: paramagnet ρ = I/2^N ⇒ S_αα(q) = 1 (independent of q and N)",
                    ],
                    fetch_kw=(; q=q_val, beta=HIGH_T_BETA),
                )
            end
        end
    end
end
