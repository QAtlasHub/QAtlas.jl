# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_s1heisenberg1d_obc_susc_batch.jl
#
# Susceptibility T → 0 limits for the OBC S=1 Heisenberg chain at even N:
# the gapped Haldane GS at even N is a unique S_total = 0 singlet, so
# χ_αα = β · Var(S^α_total) = 0 exactly.
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "S1Heisenberg1D — Susceptibility{XX,YY,ZZ}/OBC at even N, T→0 (#381 batch)" begin
    LOW_T_BETA = 1e6
    for J in (0.5, 1.0, 2.0)
        for N in (4, 6)  # even N + S=1 ⇒ unique singlet GS (Haldane phase)
            for q in (SusceptibilityXX(), SusceptibilityYY(), SusceptibilityZZ())
                verify(
                    S1Heisenberg1D(),
                    q,
                    OBC(N);
                    route=:second_closed_form,
                    independent=0.0,
                    agree_within=1e-8,
                    refs=[
                        "S1Heisenberg1D OBC even N: unique gapped Haldane S_total=0 GS ⇒ χ_αα = β·Var(S^α_total) = 0",
                    ],
                    fetch_kw=(; J=J, beta=LOW_T_BETA),
                )
            end
        end
    end
end
