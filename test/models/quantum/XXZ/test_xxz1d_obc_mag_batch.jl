# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/XXZ/test_xxz1d_obc_mag_batch.jl
#
# Symmetry-zero verification cards for the OBC XXZ1D chain: the total
# magnetisations <S^α>/N (α ∈ {x, y, z}) vanish in the thermal ensemble
# at any temperature and parameters.
#   * <S^z> = 0 by U(1) conservation of total S^z + half-filling sector.
#   * <S^x>, <S^y> = 0 by U(1) z-axis rotation symmetry (thermal trace
#     averages over the rotation, killing the x/y components).
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "XXZ1D — Magnetization{X,Y,Z}/OBC = 0 (U(1) + Z2) (#381 batch)" begin
    for (J, Δ) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (2.0, -0.5))
        for N in (4, 6, 8)
            for β in (1.0, 10.0)
                for q in (MagnetizationX(), MagnetizationY(), MagnetizationZ())
                    verify(
                        XXZ1D(),
                        q,
                        OBC(N);
                        route=:second_closed_form,
                        independent=0.0,
                        agree_within=1e-10,
                        refs=["XXZ1D U(1) z-rotation symmetry + total-S^z conservation ⇒ <S^α>_β = 0 (α ∈ {x,y,z}) in the symmetric thermal ensemble"],
                        fetch_kw=(; J=J, Δ=Δ, beta=β),
                    )
                end
            end
        end
    end
end
