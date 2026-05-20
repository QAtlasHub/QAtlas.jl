# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_heisenberg1d_obc_massgap_batch.jl
#
# At odd N the spin-1/2 OBC Heisenberg chain has a unique S_total = 1/2
# doublet ground state with two degenerate m_z = ±1/2 levels, so the
# spectral gap above the GS sector vanishes identically: Δ = 0.
# (Even N gives a finite singlet-triplet gap depending on N, not covered
# here.)  Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Heisenberg1D — MassGap/OBC at odd N = 0 (doublet GS) (#381 batch)" begin
    for J in (0.5, 1.0, 2.0)
        for N in (3, 5, 7)  # odd ⇒ S_total = 1/2 doublet ⇒ Δ = 0
            # evals[2]-evals[1] on the degenerate doublet gives machine-precision zero
            verify(
                Heisenberg1D(),
                MassGap(),
                OBC(N);
                route=:second_closed_form,
                independent=0.0,
                agree_within=1e-10,
                refs=["Heisenberg1D OBC odd N: S_total=1/2 doublet GS (two degenerate m_z = ±1/2 levels) ⇒ Δ = 0 exactly"],
                fetch_kw=(; J=J),
            )
        end
    end
end
