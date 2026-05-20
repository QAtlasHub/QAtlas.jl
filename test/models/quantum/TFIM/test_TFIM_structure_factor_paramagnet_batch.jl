# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_structure_factor_paramagnet_batch.jl
#
# At T → ∞ the TFIM thermal state is the maximally mixed state
# ρ = I / 2^N, which gives ⟨σ^α_i σ^α_j⟩ = δ_{ij} (uncorrelated) and
# therefore S_αα(q) = Σ_r e^{iqr} ⟨σ^α_0 σ^α_r⟩ → 1 for all q, all axes,
# any BC. Pure verify(); branches off main. Refs #381.
#
# TODO(#438): This card uses the σ-matrix convention (⟨σ^α⟩² = 1 ⇒
# independent = 1.0). Under the project-wide spin-S convention policy
# being rolled out in PR #438 (docs/src/conventions.md), TFIM observable
# returns must be reported in S = σ/2 units, which would change the
# expected paramagnet value to ⟨S^α_i S^α_j⟩|_{i=j} = 1/4. Migrate this
# file (independent = 0.25, tighten agree_within to 1e-9) once #438
# merges; until then the σ convention is internally consistent with the
# current TFIM structure-factor implementation.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — {XX,YY,ZZ}StructureFactor at T→∞ paramagnet = 1 (#381 batch)" begin
    HIGH_T_BETA = 1e-3

    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (0.5, 2.0))
        for q in (XXStructureFactor(), YYStructureFactor(), ZZStructureFactor())
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
end
