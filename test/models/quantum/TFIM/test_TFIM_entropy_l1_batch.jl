# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_entropy_l1_batch.jl
#
# Single-site (ℓ=1) entanglement entropy of the OBC TFIM at two trivial
# parameter limits:
#   * J = 0 (pure transverse field, GS = product state |+⟩^N at T → 0):
#       ρ₁ = |+⟩⟨+| (pure) ⇒ S_vN(ℓ=1) = S_α(ℓ=1) = 0
#   * h = 0 (pure Ising, Z2 symmetry σ_x ↔ σ_x not, σ_z ↔ −σ_z):
#       ρ₁ = I/2 (Z2-symmetric thermal trace) ⇒ S(ℓ=1) = log 2 for any T
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — VonNeumannEntropy + RenyiEntropy/OBC at ℓ=1 trivial limits (#381 batch)" begin
    # J = 0 at T → 0: pure product GS ⇒ S(ℓ=1) = 0
    for h in (0.5, 1.0, 2.0)
        for N in (4, 6, 8)
            verify(
                TFIM(; J=0.0, h=h),
                VonNeumannEntropy(),
                OBC(N);
                route=:second_closed_form,
                independent=0.0,
                agree_within=1e-10,
                refs=["TFIM J=0 T→0: GS is pure product state |+⟩^N ⇒ ρ₁ pure ⇒ S_vN(ℓ=1) = 0"],
                fetch_kw=(; ℓ=1, beta=1e6),
            )
            for α in (2, 3)
                verify(
                    TFIM(; J=0.0, h=h),
                    RenyiEntropy(α),
                    OBC(N);
                    route=:second_closed_form,
                    independent=0.0,
                    agree_within=1e-10,
                    refs=["TFIM J=0 T→0: pure GS ⇒ S_α(ℓ=1) = 0 for any Rényi index"],
                    fetch_kw=(; ℓ=1, beta=1e6),
                )
            end
        end
    end

    # h = 0 (pure Ising): Z2 symmetry σ_z → −σ_z keeps ⟨σ_z⟩ = 0, and there's
    # no x-component, so ρ₁ = I/2 ⇒ S(ℓ=1) = log 2 for any β.
    for J in (0.5, 1.0, 2.0)
        for N in (4, 6, 8)
            for β in (0.5, 10.0, 1e6)
                verify(
                    TFIM(; J=J, h=0.0),
                    VonNeumannEntropy(),
                    OBC(N);
                    route=:second_closed_form,
                    independent=log(2),
                    agree_within=1e-10,
                    refs=["TFIM h=0: Z2 symmetry σ_z → −σ_z + no x channel ⇒ ρ₁ = I/2 ⇒ S_vN(ℓ=1) = log 2 for any β"],
                    fetch_kw=(; ℓ=1, beta=β),
                )
                for α in (2, 3)
                    verify(
                        TFIM(; J=J, h=0.0),
                        RenyiEntropy(α),
                        OBC(N);
                        route=:second_closed_form,
                        independent=log(2),
                        agree_within=1e-10,
                        refs=["TFIM h=0: ρ₁ = I/2 ⇒ S_α(ℓ=1) = log 2 for any Rényi index, any β"],
                        fetch_kw=(; ℓ=1, beta=β),
                    )
                end
            end
        end
    end
end
