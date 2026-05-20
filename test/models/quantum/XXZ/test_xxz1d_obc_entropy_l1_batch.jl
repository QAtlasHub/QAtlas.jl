# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/XXZ/test_xxz1d_obc_entropy_l1_batch.jl
#
# Single-site entanglement entropy of the OBC spin-1/2 XXZ chain at
# half-filling: U(1) z-rotation symmetry + xy-Z2 reflection forces the
# single-site reduced density matrix to be ρ₁ = I/2 (no preferred axis,
# ⟨S^z⟩ = 0 in m_z=0 sector), so S_vN(ℓ=1) = S_α(ℓ=1) = log 2 for ALL
# J, Δ, N, β.
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "XXZ1D — VonNeumannEntropy + RenyiEntropy/OBC at ℓ=1 (#381 batch)" begin
    for (J, Δ) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (2.0, -0.5))
        for N in (3, 4, 5, 6)
            for β in (0.5, 10.0, 1e6)
                verify(
                    XXZ1D(),
                    VonNeumannEntropy(),
                    OBC(N);
                    route=:second_closed_form,
                    independent=log(2),
                    agree_within=1e-10,
                    refs=["XXZ1D U(1)+xy-Z2 symmetry: ρ₁ = I/2 ⇒ S_vN(ℓ=1) = log 2 for all J, Δ, N, β"],
                    fetch_kw=(; J=J, Δ=Δ, ℓ=1, beta=β),
                )
                for α in (2, 3)
                    verify(
                        XXZ1D(),
                        RenyiEntropy(α),
                        OBC(N);
                        route=:second_closed_form,
                        independent=log(2),
                        agree_within=1e-10,
                        refs=["XXZ1D U(1)+xy-Z2: maximally mixed ρ₁ ⇒ S_α(ℓ=1) = log 2 for all Rényi index"],
                        fetch_kw=(; J=J, Δ=Δ, ℓ=1, beta=β),
                    )
                end
            end
        end
    end
end
