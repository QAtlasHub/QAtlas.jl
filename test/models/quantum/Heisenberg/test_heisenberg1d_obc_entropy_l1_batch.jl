# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_heisenberg1d_obc_entropy_l1_batch.jl
#
# Single-site entanglement entropy of the OBC spin-1/2 Heisenberg chain
# is parameter-/temperature-independent: by SU(2) symmetry the reduced
# density matrix of any single site is ρ₁ = I/2 (full SU(2) average over
# orientations), so S_vN(ℓ=1) = S_α(ℓ=1) = log 2 for ALL J, N, β.
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Heisenberg1D — VonNeumannEntropy + RenyiEntropy/OBC at ℓ=1 (#381 batch)" begin
    for J in (0.5, 1.0, 2.0)
        for N in (3, 4, 5, 6)
            for β in (0.5, 10.0, 1e6)
                verify(
                    Heisenberg1D(),
                    VonNeumannEntropy(),
                    OBC(N);
                    route=:second_closed_form,
                    independent=log(2),
                    agree_within=1e-10,
                    refs=["Heisenberg1D SU(2) symmetry: ρ₁ = I/2 ⇒ S_vN(ℓ=1) = log 2 for all J, N, β"],
                    fetch_kw=(; J=J, ℓ=1, beta=β),
                )
                for α in (2, 3)
                    verify(
                        Heisenberg1D(),
                        RenyiEntropy(α),
                        OBC(N);
                        route=:second_closed_form,
                        independent=log(2),
                        agree_within=1e-10,
                        refs=["Heisenberg1D SU(2): ρ₁ = I/2 (maximally mixed) ⇒ S_α(ℓ=1) = log 2 for all Rényi index α"],
                        fetch_kw=(; J=J, ℓ=1, beta=β),
                    )
                end
            end
        end
    end
end
