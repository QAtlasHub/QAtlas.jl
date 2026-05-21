# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_fidelity_susc_J0_batch.jl
#
# At J = 0 the TFIM Hamiltonian H = -h Σ σ_x has GS = |+⟩^N (trivial
# paramagnet) independent of h — only the energy depends on h, not the
# eigenstate.  The fidelity susceptibility χ_F here is computed with
# respect to the transverse field h (∂_h-derivative); since ⟨ψ(h)|ψ(h+dh)⟩ = 1
# identically, χ_F(∂_h) = 0 exactly.  Note: this is NOT the same as ∂_J χ_F,
# which is non-zero at J = 0.
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
                refs=[
                    "χ_F via ∂_h (transverse field derivative): at J=0 the GS is the trivial paramagnet |+⟩^N independent of h, so ⟨ψ(h)|ψ(h+dh)⟩ = 1 ⇒ χ_F(∂_h) = 0 exactly. (∂_J derivative would be non-zero — convention here is ∂_h.)",
                ],
            )
        end
    end
end
