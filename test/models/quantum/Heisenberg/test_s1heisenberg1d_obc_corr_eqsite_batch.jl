# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_s1heisenberg1d_obc_corr_eqsite_batch.jl
#
# Equal-site spin-spin correlator of the OBC spin-1 Heisenberg chain:
#   ⟨S^α_i S^α_i⟩_β = ⟨(S^α)²⟩_β
# By SU(2) symmetry the single-site reduced state is ρ_i = I/3, so
# ⟨(S^α)²⟩ = (1/3) Tr(S_α²) = S(S+1)/3 = 2/3 for S=1, independent of
# the axis α, J, N, β.
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
                        refs=["S1Heisenberg1D SU(2): ρ_i = I/3 ⇒ ⟨(S^α)²⟩ = S(S+1)/3 = 2/3, axis- and β-independent"],
                        fetch_kw=(; J=J, beta=β, i=1, j=1),
                    )
                end
            end
        end
    end
end
