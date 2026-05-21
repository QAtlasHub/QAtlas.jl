# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_heisenberg1d_obc_susc_batch.jl
#
# Susceptibility T → 0 limits for the OBC Heisenberg chain at even N:
# the ground state is a unique S_total = 0 singlet, so the linear
# response χ_αα = β · Var(S^α_total) vanishes identically (Var = 0 on
# a singlet). At odd N the GS is a doublet (S_total = 1/2) and χ
# diverges as ~β/N — handled separately.
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Heisenberg1D — Susceptibility{XX,YY,ZZ}/OBC at even N, T→0 (singlet GS) (#381 batch)" begin
    LOW_T_BETA = 1e6
    for J in (0.5, 1.0, 2.0)
        for N in (4, 6, 8)  # even N ⇒ singlet GS
            for q in (SusceptibilityXX(), SusceptibilityYY(), SusceptibilityZZ())
                verify(
                    Heisenberg1D(),
                    q,
                    OBC(N);
                    route=:second_closed_form,
                    independent=0.0,
                    agree_within=1e-8,
                    refs=[
                        "Heisenberg1D OBC even N: unique S_total=0 singlet GS ⇒ χ_αα = β·Var(S^α_total) = 0",
                    ],
                    fetch_kw=(; J=J, beta=LOW_T_BETA),
                )
            end
        end
    end
end
