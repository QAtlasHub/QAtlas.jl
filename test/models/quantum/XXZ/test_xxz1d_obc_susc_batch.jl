# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/XXZ/test_xxz1d_obc_susc_batch.jl
#
# T → 0 susceptibility χ_zz of the OBC XXZ chain:
#   * Even N: unique m_z = 0 singlet GS ⇒ Var(S^z_total) = 0 ⇒ χ_zz = 0
#   * Odd N : m_z = ±1/2 doublet GS ⇒ Var(σ^z_total) = 1 ⇒
#             χ_zz = β · Var / N = β / N   (σ-matrix normalisation)
# Both branches hold for all Δ (by U(1) S^z conservation).
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "XXZ1D — SusceptibilityZZ/OBC T→0 (#381 batch)" begin
    BETA = 1e6
    for (J, Δ) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0))
        # Even N: χ_zz = 0
        for N in (4, 6, 8)
            verify(
                XXZ1D(),
                SusceptibilityZZ(),
                OBC(N);
                route=:second_closed_form,
                independent=0.0,
                agree_within=1e-9,
                refs=["XXZ1D OBC even N: singlet GS m_z=0 ⇒ Var(S^z_total)=0 ⇒ χ_zz=0"],
                fetch_kw=(; J=J, Δ=Δ, beta=BETA),
            )
        end
        # Odd N: χ_zz = β/N (σ-matrix convention)
        for N in (3, 5, 7)
            verify(
                XXZ1D(),
                SusceptibilityZZ(),
                OBC(N);
                route=:second_closed_form,
                independent=BETA / N,
                agree_within=1e-6,  # odd-N branch: independent=β/N=1e5..3e5; relative ~3e-15 floor is borderline
                refs=["XXZ1D OBC odd N: m_z=±1/2 doublet GS ⇒ χ_zz = β·Var(σ^z_total)/N = β/N"],
                fetch_kw=(; J=J, Δ=Δ, beta=BETA),
            )
        end
    end
end
