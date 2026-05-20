# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_fidelity_susc_J0_batch.jl
#
# At J = 0 the TFIM Hamiltonian H = -h Σ σ_x has GS = |+⟩^N independent of h
# (only the energy depends on h, not the eigenstate).  Changing h is a pure
# rescaling that does not rotate the GS, so the fidelity susceptibility
# χ_F = |⟨∂_h ψ_0|∂_h ψ_0⟩_⊥|² vanishes identically.
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — FidelitySusceptibility/OBC at J=0 = 0 (#381 batch)" begin
    for h in (0.5, 1.0, 2.0)
        for N in (4, 8, 12)
            verify(
                TFIM(; J=0.0, h=h),
                FidelitySusceptibility(),
                OBC(N);
                route=:second_closed_form,
                independent=0.0,
                agree_within=1e-12,
                refs=["TFIM J=0: GS = |+⟩^N is h-independent ⇒ fidelity susceptibility χ_F = 0 exactly"],
            )
        end
    end
end
