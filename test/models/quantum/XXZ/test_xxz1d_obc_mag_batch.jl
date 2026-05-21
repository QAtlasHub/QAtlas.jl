# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/XXZ/test_xxz1d_obc_mag_batch.jl
#
# Symmetry-zero verification cards for the OBC XXZ1D chain: the total
# magnetisations <S^α>/N (α ∈ {x, y, z}) vanish in the thermal ensemble
# at any temperature and parameters.
#   * <S^z> = 0 by Z₂ spin-flip symmetry K = ∏_i σ^x_i:
#         K H_XXZ K† = H_XXZ (each σ^z σ^z bond flips sign twice),
#         K S^z_tot K† = -S^z_tot ⇒ Tr(S^z_tot e^{-βH}) = 0.
#     (U(1) total-S^z conservation alone is NOT sufficient — e.g. H+h·S^z_tot
#      still has U(1) but breaks <S^z> = 0.)
#   * <S^x>, <S^y> = 0 by U(1) z-axis rotation symmetry: R_z(π) S^{x,y} R_z(π)†
#     = -S^{x,y} and [R_z, H_XXZ] = 0 ⇒ trace vanishes.
# Symmetry is BC-independent; OBC chosen because that is the registered fetch hub.
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "XXZ1D — Magnetization{X,Y,Z}/OBC = 0 (Z₂ + U(1)) (#381 batch)" begin
    for (J, Δ) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (2.0, -0.5))
        for N in (4, 6, 8)
            for β in (1.0, 10.0)
                for q in (MagnetizationX(), MagnetizationY(), MagnetizationZ())
                    verify(
                        XXZ1D(; J=J, Δ=Δ),
                        q,
                        OBC(N);
                        route=:limiting_case,
                        independent=0.0,
                        agree_within=1e-10,
                        refs=[
                            "XXZ1D OBC: Z₂ spin-flip (∏σˣ) ⇒ ⟨Sᶻ⟩ = 0; U(1) z-rotation ⇒ ⟨Sˣ⟩ = ⟨Sʸ⟩ = 0 (no field)",
                        ],
                        fetch_kw=(; beta=β),
                    )
                end
            end
        end
    end
end
