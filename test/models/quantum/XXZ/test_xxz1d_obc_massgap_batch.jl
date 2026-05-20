# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/XXZ/test_xxz1d_obc_massgap_batch.jl
#
# At odd N the spin-1/2 OBC XXZ chain has a doublet ground state in
# the m_z = ±1/2 sectors (by total-S^z parity + Z2 inversion of the xy
# plane), so the spectral gap above the GS sector vanishes: Δ = 0.
# Holds for all Δ at odd N.
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "XXZ1D — MassGap/OBC at odd N = 0 (doublet GS) (#381 batch)" begin
    for (J, Δ) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (2.0, -0.5))
        for N in (3, 5, 7)
            verify(
                XXZ1D(),
                MassGap(),
                OBC(N);
                route=:second_closed_form,
                independent=0.0,
                agree_within=1e-10,
                refs=["XXZ1D OBC odd N: m_z = ±1/2 doublet GS by U(1) S^z conservation + xy-Z2 symmetry ⇒ Δ = 0"],
                fetch_kw=(; J=J, Δ=Δ),
            )
        end
    end
end
