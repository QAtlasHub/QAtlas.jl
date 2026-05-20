# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_s1heisenberg1d_obc_corr_eqsite_batch.jl
#
# Equal-site spin-spin correlator of the OBC spin-1 Heisenberg chain:
#   ⟨S^α_i S^α_i⟩_β = ⟨(S^α)²⟩_β
# By SU(2) covariance of the thermal ensemble + spin-1 Schur's lemma:
# the rank-2 tensor T^{αβ}_i ≡ ⟨S^α_i S^β_i⟩_β commutes with global SU(2),
# so T^{αβ}_i ∝ δ^{αβ} with trace Σ_α T^{αα}_i = ⟨S_i·S_i⟩ = S(S+1) = 2,
# hence ⟨S^α_i S^α_i⟩ = ⟨(S^α)²⟩ = (1/3)·S(S+1) = 2/3 for S=1, at every
# site i, axis α, J, N, β. (Note: this does NOT require ρ_i = I/3 — that
# stronger statement fails at finite-N OBC edge sites in the Haldane regime.)
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "S1Heisenberg1D — {XX,YY,ZZ}Correlation/OBC at i=j (SU(2)) (#381 batch)" begin
    for J in (0.5, 1.0, 2.0)
        for N in (3, 4, 5, 6)
            for β in (0.5, 10.0, 1e6)
                for q in (XXCorrelation(mode=:static), YYCorrelation(mode=:static), ZZCorrelation(mode=:static))
                    verify(
                        S1Heisenberg1D(),
                        q,
                        OBC(N);
                        route=:second_closed_form,
                        independent=2/3,
                        agree_within=1e-8,
                        refs=["S1Heisenberg1D SU(2)-covariance + spin-1 Schur: ⟨S^α_i S^α_i⟩ = ⟨(S^α)²⟩ = (1/3)·S(S+1) = 2/3 at any i, axis α, J, N, β (no ρ_i = I/3 assumption — fails at finite-N OBC edge sites in Haldane regime)"],
                        fetch_kw=(; J=J, beta=β, i=1, j=1),
                    )
                end
            end
        end
    end
end
