# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/misc/test_hubbard1d_luttinger_U0_batch.jl
#
# Hubbard1D Luttinger parameter at U = 0 (free spinful fermion at half
# filling): K = 1 exactly.  Off U = 0 the model requires the Voit 1995
# Lieb–Wu integrals (Phase 2, not closed-form).
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Hubbard1D — LuttingerParameter/Infinite at U=0 = 1 (#381 batch)" begin
    for t in (0.5, 1.0, 2.0)
        verify(
            Hubbard1D(; t=t, U=0.0, μ=0.0),
            LuttingerParameter(),
            Infinite();
            route=:second_closed_form,
            independent=1.0,
            agree_within=1e-14,
            refs=[
                "At U=0 the Hubbard model is a free spinful fermion at half-filling (μ = U/2 = 0) ⇒ both Luttinger parameters K_ρ = K_σ = 1 (free-boson value, spin-charge separated free chains)",
            ],
        )
    end
end
